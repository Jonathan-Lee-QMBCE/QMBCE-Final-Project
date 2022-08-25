# This file contains all the functions necessary to run the Scalescape analysis 
# of the biological problem investigated in the project; relating species 
# richness of birds at sites across the contiguous United States to their
# surrounding landcover, as well as a random observer effect and the average
# summer temperature of the site.
#
# The current working directory *MUST* be set to the source files location.
# Raw datasets should be placed *AS DOWNLOADED* in their appropriate folders.
# The code works on the principle of lazy evaluation, where objects are loaded
# or created only when needed. All necessary pre-processing steps will be run
# automatically. Note that this may take a very long time and be processor
# intensive if steps such as converting the landcover raster need to be
# undertaken.

#Loading packages

library(sf)         # Widely used for vector objects such as sites and buffers.
library(utils)      # Used widely but infrequently for general functions.
library(tidyverse)  #
library(dplyr)      # Used for dataframe sorting such as filtering by value.
library(digest)     # Used to create hashes for creating unique filenames for objects such as dataframes
library(raster)     # Widely used for raster objects such as the NLCD cover data
library(parallel)   # Used to parallelise processes such as creating landscape matrices
library(scalescape) # Used to perform analysis, as well as bespoke data pre-processing
library(lme4)       # Used for construction of the GLMM local model used in the analysis

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

#This function finds or creates the landscape matrices for a given combination
#of sites, radius, year, and list of layers. where:
#
# sites  = A simple feature collection, containing the sites the matrix is
#          associated with. Scalescape will generated the matrix data for each
#          and relate them during model fitting by assuming they are in the same
#          order.
#
# masks  = A rasterbrick object containing the binary maps for each cover type.
#
# layers = A vector of strings, representing the cover types whose matrix for
#          the site-radius combo is required.
#
# radius = The radius around each site Scalescape should consider for matrix
#          construction. 3000 represents the upper limit of viable generation
#          times, from my experience.
#
# year   = The year the data are from. This is only used for saving/loading to
#          disk.
#
# nCores = The number of cores that the function should use to generate matrices
#          in parallel. Defaults to one less than the number available, to
#          avoid grinding the system to a halt.
#
#The function will return a list of the requested matrices, named by the name of
#their cover type.
#
#This function does not make use of the standarised "load or generate" model
#used for the rest of the code, since the large amount of time taken to generate
#the matrices incentives them to be generated and stored separately, but the
#matrices are only ever used elsewhere combined as a named list.

generateLandscapeMatrix <- function(sites,masks,layers,radius,year,nCores=detectCores()-1){
  
  #Unhelpfully, the process by which the binary cover masks are created results
  #in their having the generic names X1..8. Since they are guaranteed to be in 
  #alphabetical order due to the process which created them, it is simply but
  #clunky to relate the two. This section creates a key lookup for this, to
  #simplify the rest of the function.
  layerLookup <- new.env()
  groups <- c("Cropland","Forest","Grassland","UrbanHigh","UrbanLow","UrbanMedium","UrbanOpen","Wetland")
  for (i in 1:length(groups)){
    layerLookup[[groups[i]]] <- i
  }
  
  #pregenerates the folder where the matrices will be saved, if it does not
  #exist already.
  dirPath <-  paste("data/",year,"/matrices/r",radius,"_",digest::digest(sites$SID),sep="")
  recDir(dirPath,wd=getwd(),noFile=T)
  dirPath <- paste(getwd(),dirPath,sep="/")
  
  #layers are sorted alphabetically to ensure naming consistency.
  layers <- sort(layers)
  
  #iterate through each of the matrices requested. if they are not found on disk
  #they are added to the queue to be generated.
  pathQueue <- c()
  saveQueue <- c()
  maskList <- c()
  for (layer in layers){
    fPath <- paste(dirPath,"/",layer,".rda",sep="")
    if (file.exists(fPath) == F){
      print(paste("Landscape matrix for '",layer,"' does not exist for this site-radius combo.",sep=""))
      lk <- layerLookup[[layer]]
      maskList <- append(maskList,masks[[lk]])  
      saveQueue <- append(saveQueue,fPath)
    }
    pathQueue <- append(pathQueue,fPath)
  }
  
  #the queue of matrices to be generated is then worked through in parallel.
  #once all the matrices have been generated, they are saved to disk.
  if (length(saveQueue) > 0){
    print(paste("Generating the",length(saveQueue),"absent matri(x/ces)..."))
    cl <- makeCluster(nCores)
    landMats <- clusterMap(cl=cl,raster=maskList,fun=landscape_matrix,MoreArgs=list(sites=sites,max.radius=radius))
    stopCluster(cl)
  
    for (i in 1:length(saveQueue)){
      res <- landMats[[i]]
      file <- saveQueue[i]
      save(res,file=file)
    }
  }
  
  #now all the matrices are guaranteed to be stored on the disk, they are loaded
  #and added to the return list in alphabetical order.
  retList <- list()
  for (i in 1:length(pathQueue)){
    retList[[layers[i]]] <- get(load(pathQueue[i])[1])
  }
  return(retList)
}

#Creates a GLMM for a given set of sites, where
#
# sites = A simple feature collection, containing the dataset the model should
#         be fitted to.
#
# year  = The year the data are from. Used for saving to disk.
#
#The function doesn't actually return anything, but instead saves its result to
#disk, as per the "load-or-generate" system

generateLocalModel <- function(sites,year){
  print(paste("Generating local model for",length(sites$UID)," sites..."),sep="")
  localMod <- glmer(Richness ~ Temperature + (1 | ObsN),data=sites,family="poisson",nAGQ=2,control=glmerControl(optimizer="nloptwrap", optCtrl=list(algorithim="NLOPT_LN_NELDERMEAD")))
  savePath <- paste("data/",year,"/localMod_",digest::digest(sites),".rda",sep="")
  recDir(savePath,wd=getwd())
  save(localMod,file=savePath)
}

#Fits a Scalescape model for the given sites, considering the given land cover
#types, where:
#
# sites  = A simple feature collection, containing the dataset the model should
#          be fitted to.
#
# masks  = A rasterbrick object containing the binary maps for each cover type.
#
# year   = the year the data are from. Used for saving/loading to disk.
#
#The function doesn't actually return anything, but instead saves its result to
#disk, as per the "load-or-generate" system.

generateScalescapeModel <- function(sites,masks,layers,radius,outName,year,plotFits=F,localMod=NULL){
  layers <- sort(layers)
  
  matPath <- paste("data/",year,"/matrixList_",digest::digest(c(sites,layers,radius)),".rda",sep="")
  mats <- generateLandscapeMatrix(sites=sites,masks=masks,layers=layers,radius=radius,year=year)
  if (is.null(localMod)){
    locPath <- paste("data/",year,"/localMod_",digest::digest(sites),".rda",sep="")
    localMod <- loadOrGenerate(locPath,generateLocalModel,sites=sites,year=year)
  }
  
  landForm <- paste("~ . + ",paste(layers,collapse=" + "),sep="")
  print(paste("fitting model with landscape formula of '",landForm,"'...",sep=""))
  old <- Sys.time()
  fitted <- dist_weight(mod0 = localMod, landscape.vars = mats, landscape.formula = landForm, data = sites, plot.fits=plotFits,n.partition = 100) 
  fitted$runTime <- Sys.time()-old
  savePath <- paste("results/",outName,".rda",sep="")
  recDir(savePath,wd=getwd())
  save(fitted,file=savePath,precheck = T)
}

#Runs the Scalescape analysis, where a certain subset of the BBS data is used to
#fit a distance-weighted model relating land cover to species richness, also
#considering the observer as a random effect and the mean summer temperature of
#the site. The function will automatically generate any data it needs, assuming
#the data sources have been downloaded and placed in the correct folders.
#
#Running the function with the default values will result in the same final
#model described in the report, but values can be configured to
#types:
#
# year    = The year data should be taken from. Given the size of the NLCD data,
#           and the fact that there aren't sets for every year, this is likely
#           to be the decider of what year to use.
#
# radius  = The radius around sites that Scalescape should consider when
#           generating its landscape matrix. Generally, more is better, but a
#           value of 3000m should be more than adequate.
#
# every & = Given the size of the BBS dataset, it is best to use a subset. Every
# nth       and are nth used in conjuncation to select the subset of items whose
#           index in the set %% every == nth. For example, every = 50 and nth =
#           1 will select 1, 51, 101, etc. Given that each route is made of 50
#           stops, every = 50 is advisable, since it means the subset will
#           contain one stop per route
#
#The function will give a summary of the best model it found using the paring
#process detailed in the report, and return the model object. The model, as well
#as all the others considered during paring, will also be saved to disk.
#
#If you want to see the weighting curves for the model, the function 
#'plotScalescapeCurves' should be run with the model as the sole parameter. This
#should be used since Scalescapes automatic plot during fitting is slower and
#can sometimes crash the program.

runScaleScapeAnalysis <- function(year=2016, radius=3000, every=50, nth=1){
  
  #Sources the site data, generating it if necessary
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sites <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  
  #Sources the land cover data, generating it if necessary
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  
  #ensure sites have the same CRS as the cover data
  sitesT <- st_transform(sites,masks@crs)
  
  #subsetting
  sitesF <- filter(sitesT,row_number()%%every == nth)
  
  layers <- c("Cropland","Forest","Grassland","UrbanHigh","UrbanLow","UrbanMedium","UrbanOpen","Wetland")
  
  outName <- paste("Every",every,"_",nth,"nth",radius,"r_",paste(layers,collapse = "+"),sep="")
  loadPath <- paste("results/",outName,".rda",sep="")
  mod <- loadOrGenerate(loadPath,generateScalescapeModel,year=year,sites=sitesF,masks=masks,layers=layers,radius=radius,outName=outName,plotFits=F)

  print("full model fitted! paring down...")
  
  #Paring down works by a round based search, with each round containing simpler
  #and simpler models. The first round considers the eight simpler models
  #created by dropping one cover type from the full model. Any that outperform
  #it, judged by their having a lower AIC, will have their children, models
  #created by dropping one of their cover types, considered next round. Once
  #a round starts with no models in it, the model found that had the lowest
  #overall AIC is selected.
  
  #setup for first round
  mNames <- list()
  mLayers <- list()
  mLoadPaths <- list()
  aches <- c()
  for (i in 1:8){
    layerList <- append(head(layers,i-1),tail(layers,8-i))
    outName <- paste("Every",every,"_",nth,"nth",radius,"r_",paste(layerList,collapse = "+"),sep="")
    tLoadPath <- paste("results/",outName,".rda",sep="")
    print(i)
    mod <- loadOrGenerate(tLoadPath,generateScalescapeModel,year=year,sites=sitesF,masks=masks,layers=layerList,radius=radius,outName=outName)
    aches <- append(aches,mod$AIC)
  }
  
  lastAic <- mod$AIC
  newAic <- mod$AIC
  best <- mod
  testQueue <- dropOnePerms(layers)
  
  #since a model can be arrived at by dropping one term from several different
  #models, a blacklist is kept of models that are already being considered
  blackList <- new.env()
  
  while (length(testQueue) > 0){
    print(paste("Trying ",length(testQueue)," simpler models...",sep=""))
    newQueue <- list()
    for (i in 1:length(testQueue)){
      print(paste("",paste(testQueue[[i]],collapse="+")))
      layerList <- testQueue[[i]]
      outName <- paste("Every50_",radius,"r_",paste(layerList,collapse = "+"),sep="")
      tLoadPath <- paste("results/",outName,".rda",sep="")
      mod <- loadOrGenerate(tLoadPath,generateScalescapeModel,year=year,sites=sitesF,masks=masks,layers=layerList,radius=radius,outName=outName)
      #any model that performs better than last round's best will have its children considered next round
      if (mod$AIC < lastAic){
        if (newAic > mod$AIC){
          newAic <- mod$AIC
          best <- mod
        }
        newDrops <- dropOnePerms(names(mod$opt.range))
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
      lastAic <- newAic
    }
  }
  print("Best:")
  print(names(best$opt.range))
  print(lastAic)
  return (best)
}

#Takes a fitted Scalescape model and a new dataset of sites, and uses the model
#to predict their species richness/
#
# sSModel   = A Scalescape model produced by "generateScalescapeModel".
#
# newData   = A dataset of sites to make predictions on.
#
# masks     = A rasterbrick object containing the binary maps for each cover
#             type. These should be from the same year as the new data,
#
# year      = The year the new data are from.
#
# plotError = If set to true, the function will automatically plot the sites
#             across space, colour coded by the direction and magnitude of
#             Scalescape's prediction error/
#
#The function returns a copy of the newData dataframe, with addition columns
#"Prediction", holding the richness Scalescape predicted for each site, and
#"Error", the difference between the predicted and actual value.

makePredictions <- function(sSModel,newData,masks,year,plotError=F){
  radius <- max(sSModel$max.Dist)
  layers <- sort(names(sSModel$opt.range))
  newLandMatrices <- generateLandscapeMatrix(sites=newData,masks=masks,layers=layers,radius=radius,year=year)
  newWithLandscapes <- getWeighted(sSModel$opt.range,mod0=sSModel$mod0,landscape.formula = sSModel$landscape.formula,data=newData,max.Dist = sSModel$max.Dist, landscape.vars = newLandMatrices,weight.fn = sSModel$weight.fn)
  preds <- predict(sSModel$mod,newdata=newWithLandscapes)
  newWithLandscapes$Prediction <- exp(preds)
  newWithLandscapes$Error <- newWithLandscapes$Prediction-newWithLandscapes$Richness
  
  if (plotError){
    ggplot()+geom_sf(data=newWithLandscapes,aes(color=Error)) + scale_color_gradient2(low = "blue",mid="green", high = "red") + labs(title="Scalescape Prediction Error")
  }
  return(newWithLandscapes)
}

#Runs Scalescape's provided bootstrapping functionality on a fitted Scalescape
#model. This is *VERY SLOW* to the point of being unrunnable for the datasets
#used in this study.
#
# year    = The year the data the model is fitted to are from. Used for
#           constructing null model, as well as loading/saving to disk.
#
# model   = A Scalescape model produced by "generateScalescapeModel".
#
# outName = The name by which the results should be saved to disk.
#
#The function will return the results of the bootstrap analysis, as well as
#saving them to disk.

generateBootstrap <- function(year,model,outName){
  sites <- model$data
  print(paste("Performing bootstrap analysis on model with ",length(model$coef)-1," layers and ",length(sites$UID)," sites...",sep=""))
  locPath <- paste("data/",year,"/localMod_",digest::digest(sites),".rda",sep="")
  
  #Load or generate the NULL model that was used to create the Scalescape model
  localMod <- loadOrGenerate(locPath,generateLocalModel,sites=sites)
  
  bootMod <- dist_weight_boot(mod.full=model,mod.reduced=localMod)
  
  savePath <- paste("results/bootstraps/",outName,".rda",sep="")
  recDir(savePath,wd=getwd())
  save(bootMod,file=savePath,precheck = T)
  return(bootMod)
}

#The functions below are functions for plotting sites from the dataset. They are
#not necessary to understand the modelling, but can be very pretty. 

#Takes a site, fitted Scalescape mode, and the landcover masks, and creates a
#mosaic showing each cover type present around the site and how the model is
#weighting it. This is shown using "contour lines", with each marking a 0.1 drop
#in relevance.
plotSite <- function(site,model,masks){
  colEnv <- getColEnv()
  maskNs <- names(model$opt.range)
  ext <- extent(st_buffer(site,2000))
  maskNs <- c()
  for (name in names(model$opt.range)){
    if (sum(extract(masks[[name]],ext)) >0){
      maskNs <- append(maskNs,name)
    }
  }
  parTile(length(maskNs))
  for (maskN in maskNs){
    rads <- sqrt(-2*log(seq(0.9,0.1,-0.1)))*model$opt.range[[maskN]]
    raster::plot(masks[[maskN]],ext=ext,axes=F,legend=F,col=c("white",colEnv[[maskN]]))
    for (rad in rads){
      plot(st_buffer(site,rad),add=T,border=rgb(0,0,0,0.25),col=rgb(0,0,0.25,0.05))
    }
    plot(site,pch=3,col=rgb(1,0,0,1),cex=1,lwd=1,add=T)
    title(maskN,line=-1,lwd=2,col.main="red",bg="black")
  }
}

#for a given site and mask, displays the surroundings of the site in the given
#colour. Not helpful to the study, but sometimes very pretty.
getWideView <- function(site,mask,cola=rgb(1,0,0,1)){
  parTile(1)
  ext <- extent(st_buffer(site,5000))
  raster::plot(mask,ext=ext,axes=F,legend=F,col=c("white",cola))
  plot(site,pch=3,col=rgb(1,0,0,1),cex=1,lwd=1,add=T)
}

#semi-automatic generation of the above. randomly selects a number of sites num
#and saves a plot of the cover type ct around them to disk.
plotStuff <- function(num,ct){
  recDir(paste("spreads",ct,sep="/"),wd=getwd(),noFile=T)
  colEnv <- getColEnv()
  for (i in 1:num){
    num <- floor(runif(1, min=1, max=2161))
    png(paste("spreads/",ct,num,".png",sep=""),width=1216,height=1100,units="px")
    getWideView(sitesF[num,],masks[[ct]],colEnv[[ct]])
    dev.off()
  }
}
