# This file contains all the code related to the report's attempt to implement
# a distance weighted model using the INLAbru library.

# The current working directory *MUST* be set to the source files location.
# Raw datasets should be placed *AS DOWNLOADED* in their appropriate folders.
# The code works on the principle of lazy evaluation, where objects are loaded
# or created only when needed. All necessary pre-processing steps will be run
# automatically. Note that this may take a very long time and be processor
# intensive if steps such as converting the landcover raster need to be
# undertaken.

#Loading packages

library(sf)            # Widely used for vector objects such as sites and buffers.
library(utils)         # Used widely but infrequently for general functions.
library(tidyverse)     #
library(dplyr)         # Used for dataframe sorting such as filtering by value.
library(raster)        # Widely used for raster objects such as the NLCD cover data,
library(parallel)      # Used to parallelise the creation of the distance lists,
library(exactextractr) # Used to extract raster data in radius, since it is faster and better than Raster's.
library(inlabru)       # The focus of the study. These are used to construct and fit a model using Integrated
library(INLA)          # Nested Laplace Approximation.

#'helperFunctions' contains general purpose functions, such as the function for
#finding the index of the first occurrence of a value in a vector (or similar).
#Generally, examination of these to understand the code should not be required.
source(file="util/helperFunctions.R")

#'generators_data' contains the functions for generating the data that are
#common to all of the different analyses. This file represents the bulk of the
#data pre-processing steps.
source(file="util/generators_data.R")

#Single-core function for generating distance information mappings. Given a
#pointframe of sites, a radius, and a cover raster, the function will iterate
#through each point in the list, find all non-zero cells within the given radius
#of the site, and record the x and y component of their distance from the site,
#as well as the proportion of the buffer they occupy. These data are then stored
#in an environment, with the unique spatial identifier of the site as  the key.
#This environment is then returned.
singleGenerateDisInfoMap <- function(pointFrame,radius,rast){
  library(sf)
  library(exactextractr)
  library(dplyr)
  retEnv <- new.env()
  for (i in 1:length(pointFrame$SID)){
    point <- pointFrame[i,]
    buf <- st_buffer(point,radius)
    cells <- exact_extract(rast,buf,include_xy=T)[[1]]
    tval <- filter(cells,value==1)
    #In the event that there are no cells of the relevant type, the integration
    #will always evaluate to zero. The function therefore records a Boolean
    #FALSE value instead, which communicates that the integration should be
    #skipped.
    if (length(tval$value) == 0){
      retEnv[[point$SID]] <- list(dxs=F,dys=F,areas=F)#areas doesn't actually
    }
    else{
      distList <- rep(0,length(tval$x))
      xDistList <- rep(0,length(tval$x))
      yDistList <- rep(0,length(tval$x))
      #Get the absolute X and Y distance from the site for each cell.
      stcP <- st_coordinates(point)
      px <- stcP[[1]]
      py <- stcP[[2]]
      for (j in 1:length(tval$x)){
        xDistList[j] <- abs(tval$x[j]-px)
        yDistList[j] <- abs(tval$y[j]-py)
      }
      pList <- list()
      pList[['areas']] <- tval$coverage_fraction/sum(cells$coverage_fraction)
      pList[['dxs']] <- xDistList
      pList[['dys']] <- yDistList
      retEnv[[point$SID]] <- pList
    }
  }
  return(retEnv)
}

#Generalised version of the above, mostly for offering optional multi-threading.
#Arguments are identical, except and extra optional term dictating the number of
#cores to use, and output is identical.
generateDisInfoMap <- function(pointFrame,radius,rast,nCores=1){
  if (nCores==1){
    return (singleGenerateDisInfoMap(pointFrame,radius,rast))
  }
  #Dividing the dataframe between available cores.
  jobLots <- naiveBalance(length(pointFrame$SID),nCores)
  start <- 1
  frameList <- list()
  for (i in 1:length(jobLots)){
    end <- (start+jobLots[i]-1)
    frameList[[i]] <- pointFrame[c(start:end),]
    start <- start+jobLots[i]
  }
  cl <- makeCluster(nCores)
  newSims <- clusterMap(cl=cl,pointFrame=frameList,fun=singleGenerateDisInfoMap,MoreArgs=list(radius=radius,rast=rast))
  stopCluster(cl) 
  #Recombination.
  retEnv <- newSims[[1]]
  for (i in 2:nCores){
    tEnv <- newSims[[i]]
    for (name in names(tEnv)){
      retEnv[[name]] <- tEnv[[name]]
    }
  }
  return(retEnv)
}

#A wrapped for the "distlist" functions, abstracting away finding/generating the
#raster layer and saving the results to disk as standard in the
#"load-or-generate" system the rest of the projects code uses.
createCoverEnv <- function(year,radius,layer,nCores=1){
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sites <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  groups <- c("Cropland","Forest","Grassland","UrbanHigh","UrbanLow","UrbanMedium","UrbanOpen","Wetland")
  tarLayer <- masks[[findIn(layer,groups)]]
  
  sitesF <- st_transform(sites,masks@crs)
  
  digEnv <- generateDisInfoMap(sitesF,radius,tarLayer,nCores)
  
  path <- paste("data/",year,"/distLists/r",radius,"/",layer,".rda",sep="")
  recDir(path,wd = getwd())
  save(digEnv,file=path)
}

# All the code below is based heavily on code written by Andrew Seaton to help
# with the project. Although it has been changed quite substantially, it is
# still the same core structure, and more representative of his work than mine,
# so credit and thanks to Andy!

#The bespoke integration function, designed to operate off the digests created
#by the functions above. Takes in the X,Y,Area lists of cells around a point,
#weights them according to the (exponent of the) sigma provided, and returns the
#total.
customIntegrate <- function(lsigma,disInfo){
  sigma <- exp(lsigma)
  #if there are no cells within range, skip
  if (is.logical(disInfo[['dxs']]) == T){
    return(0)
  }
  tot <- 0
  dxs <- disInfo[['dxs']]
  dys <- disInfo[['dys']]
  are <- disInfo[['areas']]
  #2D Gaussian
  for (i in 1:length(dxs)){
    tot <- tot + are[i]*exp(-1/2*((dxs[i]/sigma)^2+(dys[i]/sigma)^2))
  }
  return(tot)
}

#A function interfacing between INLAbru and the integration function. The main
#purpose of this is to use the environment containing the data as a parameter,
#instead of a dataframe like INLAbru wants. Since INLAbru always wants to
#evaluate the same sites in the same order, this can be done safely. Returns a
#vector of values, representing each site in the data, weighted using the
#specified sigma.
kernInt <- function(sig,keyFrame){
  pSig <- paste(sig)
  if (is.null(doneEnv[[pSig]])==F){
    return(doneEnv[[pSig]])
  }
  acc <- c()
  for (i in 1:length(keyFrame$SID)){
    acc <- append(acc,customIntegrate(sig,distLists[[keyFrame$SID[i]]]))
  }
  doneEnv[[pSig]] <- acc
  return(acc)
}

doInlaBruSimulation <- function(){
  year <- 2016
  
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  sites <- loadOrGenerate(sPath,generateTemperaturedSiteList,year=year)
  
  mPath <- paste("data/",year,"bitmaps.gri",sep="/")
  masks <- loadOrGenerate(mPath,generateBinaryRasters,year=year,loadFunction=raster::brick)
  
  every <- 50
  sitesT <- st_transform(sites,masks@crs)
  sitesF <- filter(sitesT,row_number()%%every == 3)
  sitesF$Richness <- rep(0,length(sitesF$Richness))
  
  lPath <- "data/2016/distLists/r3000/Cropland.rda"
  distLists <- loadOrGenerate(lPath,genFun = createCoverEnv, year=2016,radius=3000,layer="Cropland")
  
  #DEFINE TEST SIGMA
  tesSig <- log(800)
  tesMul <- 20
  icpt <- 1
  
  #override the true species richness with a result derived from the known
  #values
  for (i in 1:length(sitesF$SID)){
    lDis <- distLists[[sitesF$SID[i]]]
    tval <- exp((tesMul*customIntegrate(tesSig,lDis))+icpt)
    sitesF$DRichness[i] <- tval
    sitesF$Richness[i] <- rpois(1,tval)
  }
  plot(sitesF$Richness~sitesF$DRichness)
  sitesSpatial <- as_Spatial(sitesF)
  
  #Convert dist env to DataFrame.
  distDatF <- data.frame(SID=sitesF$SID,disDat=rep(0,length(sitesF$SID)))
  
  #INLAbru often seems to try the same parameter values, so keeping a store of
  #sigmas and their resultant site weightings can speed up the process.
  doneEnv <- new.env()

  cmp = ~ Intercept(1) + logh(1) + beta_crop(1)
  
  fml = Richness ~ Intercept + beta_crop * kernInt(logh, .data.)
  
  lik = like(formula = fml,
             family = "poisson",
             data = distDatF,
             response_data = sitesSpatial)
  
  #Although structurally correct, this code will not arrive at any coherent
  #values. Tweaking the operating values of the "bru" command would likely be
  #key, but I did not have the knowledge or time to figure out how.
  
  fit = bru(lik,
            components = cmp,
            options = list(
              control.inla = list(int.strategy = "eb",
                                  strategy = "laplace",
                                  fast = F),
              bru_verbose = TRUE,
              bru_max_iter = 80))
  
  summary(fit)
  inlabru:::make_track_plots(fit)
}