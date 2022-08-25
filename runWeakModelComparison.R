# This file contains all the functions necessary to run the comparisons between 
# Scalescape and traditional buffer models.
#
# The current working directory *MUST* be set to the source files location.
# Raw datasets should be placed *AS DOWNLOADED* in their appropriate folders.
# The code works on the principle of lazy evaluation, where objects are loaded
# or created only when needed. All necessary pre-processing steps will be run
# automatically. Note that this may take a very long time and be processor
# intensive if steps such as converting the landcover raster need to be
# undertaken.

#Loading packages

library(geos)       # For use of its interpolate function, which is superior to raster's. #NOT HERE
library(sf)            # Widely used for vector objects such as sites and buffers.
library(utils)         # Used widely but infrequently for general functions.
library(tidyverse)     #
library(dplyr)         # Used for dataframe sorting such as filtering by value.
library(digest)        # Used to create hashes for creating unique filenames for objects such as dataframes
library(raster)        # Widely used for raster objects such as the NLCD cover data
library(parallel)      # Used to parallelise processes such as creating buffered data for sites
library(scalescape)    # Used to perform analysis, as well as bespoke data pre-processing
library(lme4)          # Used for construction of the GLMM local model used in the analysis
library(exactextractr) # Used to extract raster data in radius, since it is faster and better than Raster's.

#'helperFunctions' contains general purpose functions, such as the function for
#finding the index of the first occurrence of a value in a vector (or similar).
#Generally, examination of these to understand the code should not be required.
source(file="util/helperFunctions.R")

#'generators_data' contains the functions for generating the data that are
#common to all of the different analyses. This file represents the bulk of the
#data pre-processing steps.
source(file="util/generators_data.R")

#'scalescapeDerived' contains functions that were either directly copied from
#the 'Scalescape' package (https://github.com/benjaminiuliano/scalescape), or
#were based so heavily on its code as to be derivative. The package is licensed
#with a GPL 3 license, and is therefore permissible to edit, but such code was
#quarantined from non-derivative code for clarity.
source(file="util/scalescapeDerived.R")

#'runScalescapeAnalysis' contains all the functionality specific to Scalescape
#model construction. This is made use of in the weak model comparisons.
source(file="runScalescapeAnalysis.R")

#'runInlaBruAnalysis' contains all the functionality specific to INLAbru model
#construction. The simulated dataset created for the INLAbru modelling test is
#reused here for the simulated comparison, as is the integrate function.
source(file="runInlaBruAnalysis.R")

#Returns buffer data for the zero-radius edge case. For a given simple features
#object "sites", it will iterate through each of the masks "masks" until it
#finds the one with a non-zero value at the site, at which point it will set the
#buffer value for that cover type to 1, and all the others to zero. Radius is
#totally ignored, since it is only there to ensure that the function takes in
#the same inputs as "singleGenerateBufferedSites", so that the two can be
#interchanged easily. Returns the sites with their buffer values updated.

getCoverAtPoints <- function(sites,radius,masks){
  mapEnv <- new.env()
  for (name in names(masks)){
    mapEnv[[name]] <- paste("Buffer",name,sep="_")
  }
  for (i in 1:length(sites$SID)){
    for (cover in names(masks)){
      if (extract(masks[[cover]],sites[i,]) ==1){
        sites[[mapEnv[[cover]]]][i] <- 1
        break
      }
    }
  }
  return(sites)
}

#Takes a simple feature list "sites" and adds values to each site representing
#what proportion of the area within the radius "radius" is of that type by
#sampling a rasterbrick of cover types "masks" within that radius. Returns the
#sites with their buffer values updated.

singleGenerateBufferedSites <- function(sites,radius,masks){
  for (i in 1:length(sites$SID)){
    #exact_extract is used, since the 30mx30m resolution is substantial enough
    #at small radii to lead to unpredictable results if buffer is approximated
    #by the raster cells.
    covs <- exact_extract(masks,st_buffer(sites[i,]$Position,radius))[[1]]
    cTList <- list()
    tot <- 0
    for (cov in names(masks)){
      cSum <- sum(covs[[cov]]*covs$coverage_fraction)
      tot <- tot + cSum
      cTList[[cov]] <- cSum
    }
    for (cov in names(masks)){
      cn <- paste("Buffer",cov,sep="_")
      sites[[cn]][i] <- cTList[[cov]]/tot
    }
  }
  return(sites)
}

#generalised version of the buffering process, abstracting away the difficulties
#that come from multithreading and the possibility of a zero radius. controlled
#by arguments:
#
# year   = The year to generate buffered versions of the BBS sites for.
#
# radius = The radius of the buffer.
#
# nCores = The number of cores to use. Multithreading significantly speeds up
#          the aggregation process for large radii, but might actually be slower
#          for smaller ones.
#
#In keeping with the "load-or-generate" system, the function returns nothing,
#instead saving its output to disk.

generateBufferedSites <- function(year,radius,nCores=1){
  print(paste("Attaching ",radius,"m buffer data to sites",sep=""))
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  names(masks) <- c("Cropland","Forest","Grassland","UrbanHigh","UrbanLow","UrbanMedium","UrbanOpen","Wetland")
  
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sites <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  
  #create new columns for buffer data
  for (ct in names(masks)){
    sites[[paste("Buffer",ct,sep="_")]] <- rep(0,length(sites$SID)) 
  }
  
  sites <- st_transform(sites,masks@crs)
  
  #if the radius is zero, a different function needs to be used.
  gBFun <- singleGenerateBufferedSites
  if (radius == 0){
    gBFun <- getCoverAtPoints
  }
  
  #simple singlethreaded version
  if (nCores==1){
    sites <- gBFun(sites,radius,masks)
  }
  #more complicated multithreaded version
  else{
    #divide dataframe up between the cores
    jobLots <- naiveBalance(length(sites$SID),nCores)
    start <- 1
    frameList <- list()
    for (i in 1:length(jobLots)){
      end <- (start+jobLots[i]-1)
      frameList[[i]] <- sites[c(start:end),]
      start <- start+jobLots[i]
    }
    cl <- makeCluster(nCores)
    clusterEvalQ(cl, library("exactextractr"))
    clusterEvalQ(cl, library("raster"))
    clusterEvalQ(cl, library("sf"))
    sWBs <- clusterMap(cl=cl,sites=frameList,fun=gBFun,MoreArgs=list(radius=radius,masks=masks))
    stopCluster(cl) 
    
    #recombine outputs into a single dataframe again
    sites <- sWBs[[1]]
    for (i in 2:nCores){
      sites <- rbind(sites,sWBs[[i]])
    }
  }
  sPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",radius,".rda",sep="")
  recDir(sPath,wd=getwd())
  save(sites,file=sPath,precheck = T)
}

#This function will fit the zero-radius version of the buffer model, where only
#land cover at each sites exact location is considered. As standard, the model
#is a GLMM which includes a value for summer temperature, as well as a random
#effect for the surveyor. Parameters are:
#
# year  = The year the BBS and land cover data should be taken from.
#
# every = How the sites should be subsetted. Sites where their index %% every =
# & nth   nth will be used for fitting. For example, every = 50 and nth = 1
#         would result in sites 1, 51, 101, etc. Since every 50 sites represent
#         stops along 1 route, keeping every at 50 means that one stop from each
#         route will be used for fitting.
#
#As standard, the function returns no value, instead saving its output to disk.

runPointModelAnalysis <- function(year,every=50,nth=1){
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",radius,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=radius,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == nth)
  
  mod <- glmer(Richness ~ Temperature + CoverAtPoint + (1 | ObsN),data=sitesSub,family="poisson",nAGQ=2,control=glmerControl(optimizer="nloptwrap", optCtrl=list(algorithim="NLOPT_LN_NELDERMEAD")))
  sPath <- paste("results/weakModels/pointMod_",year,"_r",radius,".rda",sep="")
  recDir(sPath,wd=getwd())
  save(mod,file=sPath,precheck = T)
}

#Takes a buffered model "mod" and the sites it was trained on "sites", and pares
#down the model using the same method used for Scalescape models. Returns the
#best model found by the paring process, judged by AIC.

pareDownBufferModel <- function(sites,mod){
  
  #Paring down works by a round based search, with each round containing simpler
  #and simpler models. The first round considers the eight simpler models
  #created by dropping one cover type from the full model. Any that outperform
  #it, judged by their having a lower AIC, will have their children, models
  #created by dropping one of their cover types, considered next round. Once
  #a round starts with no models in it, the model found that had the lowest
  #overall AIC is selected.
  
  #First round setup.
  bufs <- c("Buffer_Cropland","Buffer_Forest","Buffer_Grassland","Buffer_UrbanHigh","Buffer_UrbanLow","Buffer_UrbanMedium","Buffer_UrbanOpen","Buffer_Wetland")
  best <- mod
  newAIC <- extractAIC(mod)[2]
  lastAIC <- newAIC
  testQueue <- dropOnePerms(bufs)
  
  #since a model can be arrived at by dropping one term from several different
  #models, a blacklist is kept of models that are already being considered
  blackList <- new.env()
  
  while (length(testQueue) > 0){
    print(paste("Trying ",length(testQueue)," simpler models...",sep=""))
    newQueue <- list()
    for (i in 1:length(testQueue)){
      bufList <- testQueue[[i]]
      form <- paste("Richness ~ Temperature",paste(bufList,collapse=" + "),"(1 | ObsN)",sep=" + ")
      mod <- glmer(form,data=sites,family="poisson",nAGQ=2,control=glmerControl(optimizer="nloptwrap", optCtrl=list(algorithim="NLOPT_LN_NELDERMEAD")))
      modAIC <- extractAIC(mod)[2]
      if (modAIC < lastAIC){
        if (newAIC > modAIC){
          newAIC <- modAIC
          best <- mod
        }
        newDrops <- dropOnePerms(bufList)
        for (j in 1:length(newDrops)){
          sn <- paste(sort(newDrops[[j]]),collapse="+")
          if (is.null(blackList[[sn]])){
            blackList[[sn]] <- T
            newQueue <- append(newQueue,list(sort(newDrops[[j]])))
          }
        }
      }
    }
    testQueue <- newQueue
    if (length(testQueue)>0){
      lastAIC <- newAIC
    }
  }
  return(best)
}

#This function will fit normal version of the buffer model, relating land cover
#in a given radius around each site to its species richness. As standard, the
#model is a GLMM which includes a value for summer temperature, as well as a
#random effect for the surveyor. Parameters are:
#
# year   = The year the BBS and land cover data should be taken from.
#
# radius = The radius around each site that should be considered,
#
# every  = How the sites should be subsetted. Sites where their index %% every =
# & nth    nth will be used for fitting. For example, every = 50 and nth = 1
#          would result in sites 1, 51, 101, etc. Since every 50 sites represent
#          stops along 1 route, keeping every at 50 means that one stop from
#          each route will be used for fitting.
#
#As standard, the function returns no value, instead saving its output to disk.

runBufferModelAnalysis <- function(year,radius=1000,every=50,nth=1){
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",radius,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=radius,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == nth)
  
  bufs <- c("Buffer_Cropland","Buffer_Forest","Buffer_Grassland","Buffer_UrbanHigh","Buffer_UrbanLow","Buffer_UrbanMedium","Buffer_UrbanOpen","Buffer_Wetland")
  form <- paste("Richness ~ Temperature",paste(bufs,collapse=" + "),"(1 | ObsN)",sep=" + ")
  mod <- glmer(form,data=sitesSub,family="poisson",nAGQ=2,control=glmerControl(optimizer="nloptwrap", optCtrl=list(algorithim="NLOPT_LN_NELDERMEAD")))
  mod <- pareDownBufferModel(sites=sitesSub,mod=mod)
  
  sPath <- paste("results/weakModels/bufferMod_",year,"_e",every,"_n",nth,"_r",r,".rda",sep="")
  recDir(sPath)
  save(mod,file=sPath,precheck = T)
}

#This function will fit the null model used for all models in the study; A GLMM
#relating species richness to summer temperature and a random effect for the
#surveyor. Parameters are:
#
# year   = The year the BBS and land cover data should be taken from.
#
# every  = How the sites should be subsetted. Sites where their index %% every =
# & nth    nth will be used for fitting. For example, every = 50 and nth = 1
#          would result in sites 1, 51, 101, etc. Since every 50 sites represent
#          stops along 1 route, keeping every at 50 means that one stop from
#          each route will be used for fitting.
#
#As standard, the function returns no value, instead saving its output to disk.

generateNullModel <- function(year, every=50, nth=1){
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",0,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=0,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == nth)
  form <- "Richness ~ Temperature + (1 | ObsN)"
  mod <- glmer(form,data=sitesSub,family="poisson",nAGQ=2,control=glmerControl(optimizer="nloptwrap", optCtrl=list(algorithim="NLOPT_LN_NELDERMEAD")))
  sPath <- paste("results/weakModels/nullMod_",year,".rda",sep="")
  recDir(sPath)
  save(mod,file=sPath,precheck = T)
}

#This function iterates across a range of radii in finer and finer detail. For
#each radius, it fits and pares down a buffer model, records the R^2 value for
#the result model, and plots it onto the graph. Parameters are:
#
# year   = The year the BBS and land cover data should be taken from.
#
# every  = How the sites should be subsetted. Sites where their index %% every =
# & nth    nth will be used for fitting. For example, every = 50 and nth = 1
#          would result in sites 1, 51, 101, etc. Since every 50 sites represent
#          stops along 1 route, keeping every at 50 means that one stop from
#          each route will be used for fitting.
#
# maxR   = The maximum radius the function should fit models for.
#
# jump   = Initially, the function will start at zero radius. It will then
# &        continually fit a model for the current radius and then increment it
# minJ     by the jump value until it reaches the maximum value, at which point,
#          it will reset to zero and begin again (radii already considered
#          won't be fitted again. when the jump value is halved such that it is
#          below minJ, the process wil halt.
#
# new    = A flag used to decide whether to create a new plot or not. setting
#          new to FALSE will cause the function to plot its results to the
#          latest plot. This can be useful if you want to manually refine the
#          search. Alternatively, "compileRSquaredPlot" can be used instead to
#          find and plot the previous generated models.
#
#The function returns nothing and doesn't save its output to disk. However, the
#models fitted for each radius are saved, and can be used for future plotting or
#other functions.

liveCreateRSquaredPlot <- function(year, every=50,nth=1,maxR=3000,jump=1500,minJ=50, new=T){
  nMPath <- paste("results/weakModels/nullMod_",year,".rda",sep="")
  nullMod <- loadOrGenerate(nMPath,generateNullModel,year=year,every=every,nth=nth)
  nullVar <- getME(nullMod,"devcomp")$cmp[["dev"]]

  blackList <- new.env()
  rList <- c()
  rSList <- c()
  
  if (new){
    plot(1, type = "n", xlab = "Buffer Radius", ylab = "R^2", xlim = c(0, maxR), ylim = c(0,1))
  }
  
  while (jump > minJ){
    for (r in seq(0,maxR,jump)){
      #given the halving increment search scheme, a blacklist is needed so that
      #values aren't plotted repeatedly.
      if (is.null(blackList[[paste(r)]])){
        blackList[[paste(r)]] <- T
        nMPath <- paste("results/weakModels/bufferMod_",year,"_e",every,"_n",nth,"_r",r,".rda",sep="")
        newMod <- loadOrGenerate(nMPath,runBufferModelAnalysis,year=year,every=every,radius=r,nth=nth)
        lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",r,".rda",sep="")
        sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=r,nCores=4)
        sitesSub <- filter(sites,row_number()%%every == nth)
        sites <- sitesSub
        rS <- MuMIn::r.squaredGLMM(newMod)[2,2]
        points(r,rS,col="red")
        rList <- append(rList,r)
        rSList <- append(rSList,rS)
      }
    }
    jump <- jump/2
  }
}

#This function searches the results folder for optimal models found for each
#radius, calculates their R^2 values, and then plots them on a graph. Parameters
#are:
#
# year   = The year the BBS and land cover data should be taken from.
#
# every  = The function will only find models fitted using the same data, ergo,
# & nth    the same year, every, and nth parameters.
#
#The function returns nothing and doesn't save its output to disk.

compileRSquaredPlot <- function(year=2016,every=50,nth=1){
  plotFrame <- data.frame("R"=NULL,"rsq"=NULL,"code"=NULL)
  codeMap <- new.env()
  header <- paste("bufferMod_",year,"_e",every,"_n",nth,"_r",sep="")
  for (savedModel in list.files("results/weakModels")){
    #since models are identified by a complicated but consistent naming scheme,
    #some byzantine sting manipulation is necessary.
    if (substr(savedModel,1,nchar(header)) == header){
      radStr <- tail(strsplit(savedModel,"_")[[1]],1)
      r <- as.numeric(substr(radStr,2,nchar(radStr)-4))
      mPath <- paste("results/weakModels",savedModel,sep="/")
      #dummy values are used, since the model is already known to exist on disk
      newMod <- loadOrGenerate(mPath,runBufferModelAnalysis,year=9999,every=1,radius=1,nth=1)
      
      cl <- c()
      for (dp in names(newMod@frame)){
        if (substr(dp,1,7)=="Buffer_"){
          ct <- substr(dp,8,nchar(dp))
          if (is.null(codeMap[[ct]])){
            if (substr(ct,1,1) == "U"){
              codeMap[[ct]] <- substr(ct,6,6)
            }
            else{
              codeMap[[ct]] <- substr(ct,1,1)
            }
          }
          cl <- append(cl,codeMap[[ct]])
        }
      }
      code <- paste(sort(cl),collapse = "")
      
      lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",r,".rda",sep="")
      sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=r,nCores=4)
      sitesSub <- filter(sites,row_number()%%every == nth)
      sites <- sitesSub
      rS <- MuMIn::r.squaredGLMM(newMod)[2,2]
      
      plotFrame <- rbind(plotFrame,data.frame("R"=r,"rsq"=rS,"code"=code))
    }
  }
  ggplot(plotFrame, aes(x=R, y=rsq, colour=code)) + geom_point(size=2)+guides(color = guide_legend(title = "Model\nComposition"))+labs(title="Optimal R-Squared Value by Radius",x="Buffer Radius",y="R-Squared") + scale_x_continuous(trans='log10')
}

#This function runs a comparison between Scalescape's model, a zero-distance
#buffer model, and a buffer model with a non-zero radius. The models are used
#to predict the species richness of novel data, and the accuracy of their
#predictions are used as a show of predictive power. Parameters are:
#
# year         = The year the BBS and land cover data should be taken from.
#
# bufferRadius = The radius the non-zero buffer model should use for its buffer.
#
#The function returns nothing and doesn't save its output to disk.

compareScaleScapeWithWeakModels <- function(year,bufferRadius=31.25){
  every <- 50
  nth <- 1
  
  #The 0 radius and buffer radius models are loaded.
  pmp <- paste("results/weakModels/bufferMod_",year,"_e",every,"_n",nth,"_r",0,".rda",sep="")
  pointMod <- loadOrGenerate(pmp,runBufferModelAnalysis,year=year,radius=0)
  
  bmp <- paste("results/weakModels/bufferMod_",year,"_e",every,"_n",nth,"_r",bufferRadius,".rda",sep="")
  bufferMod <- loadOrGenerate(bmp,runBufferModelAnalysis,year=year,radius=bufferRadius)
  
  #Novel data are sourced for both the buffered models.
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",bufferRadius,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=bufferRadius,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == 25)
  
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",0,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=0,nCores=4)
  sitesSubZero <- filter(sites,row_number()%%every == 25)
  
  #The data used for Scalescape are loaded.
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sitesS <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  
  sitesST <- st_transform(sitesS,masks@crs)
  
  #Training data for Scalescape are generated, in case the model is not already
  #present on disk.
  sitesSF <- filter(sitesST,row_number()%%every == nth)
  
  #The Scalescape model is loaded or generated. Note that unlike the buffer
  #models, Scalescape models do not have a fixed name. By convention, an "every"
  #value of 50 and "nth" value of 1 are used, but this cannot be verified
  #easily. This is why the "every" and "nth" for this function are hard-coded;
  #different values can be used, but the user must be cognisant of the need to
  #make certain checks about the data used to fit the model themselves.
  
  layers <- c("Cropland","Forest","Grassland","UrbanHigh","UrbanLow","UrbanOpen","Wetland")
  sSRadius <- 3000
  outName <- paste("Every",every,"_",sSRadius,"r_",paste(layers,collapse = "+"),sep="")
  loadPath <- paste("results/",outName,".rda",sep="")
  sSMod <- loadOrGenerate(loadPath,generateScalescapeModel,year=year,sites=sitesSF,masks=masks,layers=layers,radius=sSRadius,outName=outName,plotFits=F)
  
  #Predictions are generated for the novel data using both buffer models, and
  #their levels of error are also recorded.
  predictions <- data.frame("Pred_point"=exp(predict(pointMod,sitesSubZero)),"Pred_buffer"=exp(predict(bufferMod,sitesSub)))
  predictions$Error_point <- predictions$Pred_point-sitesSub$Richness
  predictions$Error_buffer <- predictions$Pred_buffer-sitesSub$Richness
  
  #Novel data for Scalescape are loaded and the Scalescape model is used to make
  #predictions of them. Scalescape's error is also tracked.
  sitesSF <- filter(sitesST,row_number()%%every == 25)
  
  sitesSF <- makePredictions(sSMod,sitesSF,masks,year)
  predictions$Pred_scalescape <- sitesSF$Prediction
  predictions$Error_scalescape <- sitesSF$Error
  
  #The buffered models errors for each site are compared to Scalescape's. A
  #value of True means that Scalescape had less error, a value of False means it
  #didn't.
  
  predictions$Error_comparison_point_scalescape <- abs(predictions$Error_point)>abs(predictions$Error_scalescape)
  predictions$Error_comparison_buffer_scalescape <- abs(predictions$Error_buffer)>abs(predictions$Error_scalescape)
  
  summary(predictions$Error_comparison_point_scalescape)
  summary(predictions$Error_comparison_buffer_scalescape)
  
  sitesSF$ECPS <- predictions$Error_comparison_point_scalescape
  sitesSF$ECBS <- predictions$Error_comparison_buffer_scalescape
  
  #Sometimes the buffer aggregation process can create NA cover values, so this
  #is placed in as a safety measure.
  sitesSF <- na.omit(sitesSF)
  
  #Plots are produced of the sites in space, colour-coded by whether
  #Scalescape's prediction was closer to the true richness. An overall count of
  #how many sites Scalescape predicted more accurately is displayed in the key.
  totF <- sum(sitesSF$ECPS == 0)
  totT <- length(sitesSF$ECPS)-totF
  ggplot()+geom_sf(data=sitesSF,pch=4,aes(color=ECPS))+scale_color_manual(labels = c(paste("False (",totF,")",sep=""), paste("True (",totT,")",sep="")),values = c("orange", "lightskyblue"))+guides(color = guide_legend(title = "Scalescape\nBetter?"))+labs(title="Scalescape Compared to Point Model")
  totF <- sum(sitesSF$ECBS == 0)
  totT <- length(sitesSF$ECBS)-totF
  ggplot()+geom_sf(data=sitesSF,pch=4,aes(color=ECBS))+scale_color_manual(labels = c(paste("False (",totF,")",sep=""), paste("True (",totT,")",sep="")),values = c("orange", "lightskyblue"))+guides(color = guide_legend(title = "Scalescape\nBetter?"))+labs(title="Scalescape Compared to Buffer Model")
}

#This function runs a comparison between Scalescape's model, a zero-distance
#buffer model, and a buffer model with a non-zero radius using simulated data
#derived from intergrating the Cropland surrounding each site with a known sigma
#value.
#
#The function returns nothing, but saves the models fitted on the simulated data
#to disk.

doSimulatedComparison <- function(){
  year <- 2016
  every <- 50
  nth = 3
  
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sites <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  sitesT <- st_transform(sites,masks@crs)
  sitesF <- filter(sitesT,row_number()%%every == nth)
  sitesF$Richness <- rep(0,length(sitesF$Richness))
  sitesF$DRichness <- rep(0,length(sitesF$Richness))
  
  #the digest functionality developed for the INLAbru analysis is repurposed
  #here.
  lPath <- "data/2016/distLists/r3000/Cropland.rda"
  distLists <- loadOrGenerate(lPath,genFun = createCoverEnv, year=2016,radius=3000,layer="Cropland")
  
  tesSig <- log(800)
  tesMul <- 10
  icpt <- 1
  
  seed <- 806
  set.seed(seed)
  
  #Likewise, the intergration function is reused. Incidently, this is why the
  #simulation deals in the sigma value in terms of its log. This has no benefit
  #here, but is useful for the INLAbru analysis.
  for (i in 1:length(sitesF$SID)){
    lDis <- distLists[[sitesF$SID[i]]]
    tval <- exp((tesMul*customIntegrate(tesSig,lDis))+icpt)
    sitesF$Richness[i] <- rpois(1,tval)
  }
  layers <- c("Cropland")
  outName <- paste("SIMULATED_Sig",exp(tesSig),"_m",tesMul,"_i",icpt,"_Seed",seed,sep="")
  loadPath <- paste("results/",outName,".rda",sep="")
  
  #Since there are no other variables in play, the null model used by all the
  #models is an intercept-only GLM.
  localMod <- glm(Richness ~ 1,data=sitesF,family="poisson")
  
  simMod <- loadOrGenerate(loadPath,generateScalescapeModel,year=year,sites=sitesF,masks=masks,layers=layers,radius=3000,outName=outName,plotFits=F,localMod=localMod)
  
  #Generate the true curve and then plot it on top of Scalescape's estimate. The
  #function used for the curve plot corresponds to the weight of a point
  #the given distance away, but with the entire distance in one direction, for
  #example, due north of the site.
  plotScalescapeCurves(simMod)
  curve(exp(-1/2*(x/exp(tesSig))^2),from=0,to=3000,add=T,lty=2,col="#FF000088")
  
  #Since novel simulated data would be fairly similar, the models are instead
  #made to make predictions on the same dataset used for fitting.
  sitesSP <- makePredictions(simMod,sitesF,masks,year)
  predictions <- data.frame("Pred_scalescape"=sitesSP$Prediction,"Error_scalescape"=sitesSP$Error)
 
  #r3000 buffer model
  radius <- 3000
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",radius,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=radius,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == nth)
  sitesSub$Richness <- sitesF$Richness
  
  form <- Richness ~ Buffer_Cropland
  simBufMod <- glm(form,data=sitesSub,family="poisson")
  
  outname <- paste("SIMULATED_",year,"_r",radius,"_Sig",exp(tesSig),"_m",tesMul,"_i",icpt,"_Seed",seed,sep="")
  sPath <- paste("results/weakModels/",outname,".rda",sep="")
  recDir(sPath,wd = getwd())
  save(simBufMod,file=sPath,precheck = T)
  
  predictions$Pred_buffer <- exp(predict(simBufMod,sitesSub))
  predictions$Error_buffer <- predictions$Pred_buffer-sitesSub$Richness
  predictions$Error_comparison_buffer_scalescape <- abs(predictions$Error_buffer)>abs(predictions$Error_scalescape)
  
  #r0 buffer model
  lPath <- paste("data/BBS/",year,"/sites_with_AggregatedData_r",0,".rda",sep="")
  sites <- loadOrGenerate(lPath,generateBufferedSites,year=year,radius=0,nCores=4)
  sitesSub <- filter(sites,row_number()%%every == nth)
  sitesSub$Richness <- sitesF$Richness
  
  form <- Richness ~ Buffer_Cropland
  simZeroBufMod <- glm(form,data=sitesSub,family="poisson")
  
  predictions$Pred_point <- exp(predict(simZeroBufMod,sitesSub))
  predictions$Error_point <- predictions$Pred_point-sitesSub$Richness
  predictions$Error_comparison_point_scalescape <- abs(predictions$Error_point)>abs(predictions$Error_scalescape)
  outname <- paste("SIMULATED_",year,"_r0_Sig",exp(tesSig),"_m",tesMul,"_i",icpt,"_Seed",seed,sep="")
  
  sPath <- paste("results/weakModels/",outname,".rda",sep="")
  recDir(sPath,wd = getwd())
  save(simBufMod,file=sPath,precheck = T) 
  
  #The same plots are produced as in the BBS comparison: two plots of the points
  #across space, colour-coded by whether Scalescape or the respective buffer
  #model made a more accurate prediction.
  
  sitesF$ECPS <- predictions$Error_comparison_point_scalescape
  sitesF$ECBS <- predictions$Error_comparison_buffer_scalescape
  totF <- sum(sitesF$ECPS == 0)
  totT <- length(sitesF$ECPS)-totF
  ggplot()+geom_sf(data=sitesF,pch=4,aes(color=ECPS))+scale_color_manual(labels = c(paste("False (",totF,")",sep=""), paste("True (",totT,")",sep="")),values = c("orange", "lightskyblue"))+guides(color = guide_legend(title = "Scalescape\nBetter?"))+labs(title="Scalescape Compared to Point Model")
  totF <- sum(sitesF$ECBS == 0)
  totT <- length(sitesF$ECBS)-totF
  ggplot()+geom_sf(data=sitesF,pch=4,aes(color=ECBS))+scale_color_manual(labels = c(paste("False (",totF,")",sep=""), paste("True (",totT,")",sep="")),values = c("orange", "lightskyblue"))+guides(color = guide_legend(title = "Scalescape\nBetter?"))+labs(title="Scalescape Compared to Buffer Model")
}