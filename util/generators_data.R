# This file contains all the functions which generate and process the data used
# by all three analyses.
#
# The methodology used by this projects code is one of lazy evaluation. When
# something is asked for, the program checks if it is present on the disk, if
# not, it runs the function which generates it. This occurs recursively, meaning
# that, assuming the datafiles are downloaded and in the correct directories,
# any function can be run without worrying about pre-preparation. in practice,
# some of the generation functions take many hours to run, so the user should be
# cognizant of just what a functions pre-requisites are...

#Loading packages

library(geos)          # For use of its interpolate function, which is superior to raster's.
library(sf)            # Widely used for vector objects such as sites and buffers.
library(utils)         # Used widely but infrequently for general functions.
library(tidyverse)     #
library(dplyr)         # Used for dataframe sorting such as filtering by value.
library(digest)        # Used to create hashes for creating unique filenames for objects such as dataframes
library(raster)        # Widely used for raster objects such as the NLCD cover data
library(parallel)      # Used to parallelise processes such as creating buffered data for sites
library(lme4)          # Used for construction of the GLMM local model used in the analysis
library(exactextractr) # Used to extract raster data in radius, since it is faster and better than Raster's.

#For a given year, this function extracts the routes relevant to the study,
#which is to say ones that are in the contiguous United States, took place in
#the relevant year, and ran without incident. The function does not return these
#values, but instead saves them to disk.

generateRouteWhitelist <- function(year){
  print("Generating whitelist of routes from weather data")
  #It would actually be quite easy to set this up to automatically download the
  #BBS data, but it was considered a little too proactive.
  if (file.exists("data/BBS/Raw/2020Release_Nor.zip") == F){
    stop("Raw BBS data not found. Please download from 'https://www.sciencebase.gov/catalog/file/get/5ea04e9a82cefae35a129d65' and place zip file in 'data/BBS/Raw'.")
  }
  unzip("data/BBS/Raw/2020Release_Nor.zip",c("Weather.zip"),exdir="data/BBS/Raw/temp")
  unzip("data/BBS/Raw/temp/Weather.zip",c("weather.csv"),exdir="data/BBS/Raw/temp")
  weatherDat <- read.csv("data/BBS/Raw/temp/weather.csv")
  file.remove("data/BBS/Raw/temp/Weather.zip")
  file.remove("data/BBS/Raw/temp/weather.csv")
  
  #The RPID is a code expressing if anything went wrong during the survey. 101
  #means nothing did.
  qualityList <- filter(weatherDat,(CountryNum == 840) & RPID == 101 & Year == year)
  qualityList <- subset(qualityList, select = c(RouteDataID,StateNum,Route,ObsN))
  
  #Unconventional, but hopefully quite fast unique-enforcing skulduggery.
  blackList <- new.env()
  URDVec <- c()
  USNVec <- c()
  URVec <- c()
  OSNVec <- c()
  
  #Tracking duplicate runs of the same route, mostly out of sheer curiosity.
  dups <- 0
  for (i in 1:length(qualityList$Route)){
    rid <- paste(qualityList$StateNum[i],qualityList$Route[i])
    if (is.null(blackList[[rid]])){
      blackList[[rid]] <- T
      URDVec <- append(URDVec,qualityList$RouteDataID[i])
      USNVec <- append(USNVec,qualityList$StateNum[i])
      URVec <- append(URVec,qualityList$Route[i])
      OSNVec <- append(OSNVec,qualityList$ObsN[i])
    }
    else{
      dups <- dups +1
    }
  }
  UqualityList <- data.frame(RouteDataID=URDVec,StateNum=USNVec,Route=URVec,ObsN=OSNVec)
  UqualityList <- mutate(UqualityList,UID=paste(StateNum,Route,RouteDataID,sep="-"))
  
  outFile <- paste("data/BBS",year,"routes_vetted.rda",sep="/")
  if (dups > 0){
    print(paste(dups,"duplicate runs were ignored!"))
  }
  recDir(outFile,wd=getwd())
  save(UqualityList, file=outFile)
}

#This function splits each route into its 50 composite stops and combines each 
#stop with its count data, which are kept seperately in the BBS dataset. In the
#process, it also converts it to species richness. The new dataset of stops is
#saved to disk

generateStopsWithCounts <- function(year){
  print("converting whitelisted routes to stops...")
  rPath <- paste("data/BBS",year,"routes_vetted.rda",sep="/")
  routes <- loadOrGenerate(rPath,generateRouteWhitelist,year=year)
  
  print(" Preparing stops dataframe...")
  stopVec <- rep(seq(1,50),length(routes$RouteDataID))
  #Since each route contains 50 stops, the stop vector is 50 times as long.
  #"fifty-duplicated", for the record...
  FDupUID <- rep(routes$UID,each=50)
  UIDVec <- paste(FDupUID,stopVec,sep="-")
  FDupState <- rep(routes$StateNum,each=50)
  FDupObs <- rep(routes$ObsN,each=50)
  FDupPadRoutes <- rep(str_pad(routes$Route,width=3,side="left",pad="0"),each=50)
  #The BBS dataset has a standarised way of combining stop attributes to create
  #a unique(ish) code. This approach is used and extended throughout this code.
  SIDVec <- paste(FDupState,FDupPadRoutes,rep("-",length(stopVec)),stopVec,sep="")
  
  stops <- data.frame("UID"=UIDVec,"SID"=SIDVec,"Richness"=rep(NA,length(SIDVec)),"ObsN"=FDupObs)
  
  #Since the order of stops in the dataframe bears no relation to where their
  #data are in the various files the BBS data stores them in, the program
  #uses a hashed environment as a go-between to prevent lots of iteration
  #through very long lists. This also allows the program to only bother to store
  #count data from the right runs in memory.
  print(" Building count hash...")
  
  #Prime the hash with the Unique IDentifiers of each route.
  countHash <- new.env()
  for (i in 1:length(UIDVec)){
    countHash[[UIDVec[i]]] <- 0
  }
  
  #There is no plausible way this could happen, but belt and braces!
  if (file.exists("data/BBS/Raw/2020Release_Nor.zip") == F){
    stop("Raw BBS data not found. Please download and place zip file in 'data/BBS/Raw'.")
  }
  
  unzip("data/BBS/Raw/2020Release_Nor.zip",c("50-StopData.zip"),exdir="data/BBS/Raw/temp")
  
  paths <- c(NULL)
  for (i in 1:10){
    if (i == 9){
      #Bloody BBS
      fstr <- paste("50-StopData/1997ToPresent_SurveyWide/fifty",i,".zip",sep="")
    }
    else{
      fstr <- paste("50-StopData/1997ToPresent_SurveyWide/Fifty",i,".zip",sep="")  
    }
    paths <- append(paths,fstr)
  }
  
  unzip("data/BBS/Raw/temp/50-StopData.zip",paths,exdir="data/BBS/Raw/temp")
  file.remove("data/BBS/Raw/temp/50-StopData.zip")
  
  for (i in 1:10){
    print(paste(" Hashing datafile",i))
    zipPath <- paste("data/BBS/Raw/temp",paths[i],sep="/")
    fName <- paste("fifty",i,".csv",sep="")
    unzip(zipPath,c(fName),exdir="data/BBS/Raw/temp")
    file.remove(zipPath)
    path <- "data/BBS/Raw/temp"
    fPath <- paste(path,fName,sep="/")
    frame <- as.data.frame(read_csv(fPath,show_col_types = F))
    file.remove(fPath)
    frame <- filter(frame,(CountryNum == 840) & (RPID == 101) & (Year==year))
    stopCols <- colnames(frame)[8:57]
    for (j in 1:length(frame$Year)){
      uid <- paste(frame$StateNum[j],frame$Route[j],frame$RouteDataID[j],sep="-")
      #built on assumption that each species has one count per route/id/etc combo
      if (uid %in% routes$UID){
        stopNum <- 1
        for (cl in stopCols){
          hashCode <- paste(uid,stopNum,sep="-")
          chv <- countHash[[hashCode]]
          
          #One improvement explored during the study was to consider adjacent
          #stops when evaluating species richness. This was ultimately not used,
          #but two alternative aggregators have been kept below for future
          #reference.
          
          #consider 2 adjacent stops per side
          #pab <- 1*(frame[[cl]][j]>0 || (stopNum > 1 && frame[[stopCols[stopNum-1]]][j]>0) || (stopNum < 50 && frame[[stopCols[stopNum+1]]][j]>0) || (stopNum > 2 && frame[[stopCols[stopNum-2]]][j]>0) || (stopNum < 49 && frame[[stopCols[stopNum+2]]][j]>0))
          #consider adjacent stops
          #pab <- 1*(frame[[cl]][j]>0 || (stopNum > 1 && frame[[stopCols[stopNum-1]]][j]>0) || (stopNum < 50 && frame[[stopCols[stopNum+1]]][j]>0))
          #no adjacent stops
          pab <- 1*(frame[[cl]][j]>0)
          
          countHash[[hashCode]] <- chv+pab
          stopNum <- stopNum + 1
        }
      }
    }
  }
  #Once all the count data has been collected, the program iterates through the
  #stop dataframe, and reunites each stop with its count.
  print(" Adding hashed values to dataframe")
  for (i in 1:length(stops$UID)){
    stops$Richness[i] <- countHash[[stops$UID[i]]]
  }
  
  savePath <- paste("data/BBS",year,"stops_with_counts.rda",sep="/")
  recDir(savePath,wd=getwd())
  save(stops, file=savePath)
  print(" Done.")
}

#The BBS dataset contains very little spatial data about where and how the
#routes actually run. The only data of this sort is a coordinate marking the
#start of the route, which this function extracts from the dataset and saves in
#a more easy to reference format.

generateRouteStarts <- function(){
  if (file.exists("data/BBS/Raw/2020Release_Nor.zip") == F){
    stop("Raw BBS data not found. Please download and place zip file in 'data/BBS/Raw'.")
  }
  unzip("data/BBS/Raw/2020Release_Nor.zip",c("routes.zip"),exdir="data/BBS/Raw/temp")
  unzip("data/BBS/Raw/temp/routes.zip",c("routes.csv"),exdir="data/BBS/Raw/temp")
  routeDat <- read.csv("data/BBS/Raw/temp/routes.csv")
  file.remove("data/BBS/Raw/temp/routes.zip")
  file.remove("data/BBS/Raw/temp/routes.csv")
  startHash <- new.env()
  for (i in 1:length(routeDat$Route)){
    code <- paste(routeDat$StateNum[i],str_pad(routeDat$Route[i],width=3,side="left",pad="0"),sep="")
    startHash[[code]] <- st_point(c(routeDat$Longitude[i],routeDat$Latitude[i]))
  }
  savePath <- "data/Route Locations/route_starts.rda"
  recDir(savePath,wd=getwd())
  save(startHash, file=savePath)
}

#Fortunately, the USGS Patuxent Wildlife Research Center has assembled a
#database containing the spatial data for most of the BBS routes in the US. This
#function converts this into a set of points corresponding to each of a route's
#stops, and then saves them to disk.

generateStopPositions <- function(){
  sPath <- "data/Route Locations/route_starts.rda"
  routeStarts <- loadOrGenerate(sPath,generateRouteStarts)
  
  print(" Extracting spatial stop data")
  if (file.exists("data/Route Locations/Raw/Breeding Bird Survey Route Locations for Lower 48 States.zip") == F){
    stop("Raw BBS Route data not found. Please download from 'https://databasin.org/datasets/02fe0ebbb1b04111b0ba1579b89b7420/' and place zip file in 'data/Route Locations/Raw'.")
  }
  unzip("data/Route Locations/Raw/Breeding Bird Survey Route Locations for Lower 48 States.zip",exdir="data/Route Locations/Raw/temp")
  routeShapes <- read_sf("data/Route Locations/Raw/temp/data/commondata/data0/bbsrtsl020Copy.shp")
  unlink("data/Route Locations/Raw/temp",recursive = T)
  pointHash <- new.env()
  bkCount <- 0
  fails <- c()
  for (r in 1:length(routeShapes$RTENO)){
    geom <-  routeShapes$geometry[r]
    stops <- c()
    for (i in 1:50){
      #This is the main assumption made in this function. The dataset does not
      #contain positions per stop, but instead contains the entire route. In
      #theory, the stops should, as this function assumes, be evenly spaced
      #along it, but it is known this is not always the case. There is nothing
      #that can be done about this, but it is worth keeping in mind.
      stops <- append(stops,st_as_sf(geos_interpolate_normalized(geom,i-1/49)))
    }
    routeCode <- toString(routeShapes$RTENO[r])
    startPoint <- routeStarts[[routeCode]]
    if (is.null(startPoint)){
      fails <- append(fails,routeCode)
    }
    else{
      start <- st_sfc(st_cast(startPoint,"POINT"),crs="NAD83")
      bk <- F
      #Annoyingly, no indication is given in this dataset as to which end of 
      #each route is the start. This is solved by taking the start coordinate
      #specified in the BBS data, and seeing which end of the route it is closer
      #to.
      if (st_distance(stops[[50]],start) < st_distance(stops[[1]],start)){
        bkCount <- bkCount +1
        bk <- T
      }
      for (i in 1:50){
        j <- i
        if (bk){
          j <- 51-i
        }
        stopCode <- paste(routeShapes$RTENO[r],j,sep="-")
        pointHash[[stopCode]] <- stops[[i]]
      }
    }
  }
  print(paste("(",bkCount," routes were backwards!)",sep=""))
  print(paste("(",length(fails)," routes had no start point!)",sep=""))
  savePath <- "data/Route Locations/route_stops.rda"
  recDir(savePath,wd=getwd())
  #The spatial data are saved as environment, so they can be quickly reconciled
  #with their relevant stop.
  save(pointHash, file=savePath)
}

#This function simply iterates through all of the stops and attaches their
#spatial data. These new "sites" are saved to disk.

generateSiteList <- function(year){
  print("Turning stops into sites...")
  rPath <- paste("data/BBS",year,"stops_with_counts.rda",sep="/")
  sites <- loadOrGenerate(rPath,generateStopsWithCounts,year=year)
  sites$Position <- rep(NA,length(sites$UID))
  pPath <- "data/Route Locations/route_stops.rda"
  stopHash <- loadOrGenerate(pPath,generateStopPositions)
  for (i in 1:length(sites$SID)){
    sh <- stopHash[[sites$SID[i]]]
    if (is.null(sh)){
      sh <- NA
    }
    sites$Position[i] <- sh
  }
  ol <- length(sites$Position)
  #Any site where spatial data could not be extracted is stricken from the
  #dataset. In practice, this only happens on the few occasions where the UGCS
  #did not know where the route was, and so eliminates and entire routes worth
  #of stops at once.
  sites <- na.omit(sites)
  
  #Any site whose spatial data are corrupt are stricken from the study.
  #Likewise, this tends to happen on a route by route basis.
  sites <- filter(sites,st_is_empty(sites$Position)==F)
  
  print(paste("(",ol-length(sites$Position)," sites had no stop data)",sep=""))  
  sites <- st_as_sf(sites,crs="NAD83")
  savePath <- paste("data/BBS",year,"sites_no_temperature.rda",sep="/")
  recDir(savePath,wd=getwd())
  save(sites, file=savePath)
}

#This function uses the PRISM dataset to create a map of average
#summer temperatures across the United States, achieved by averaging the monthly
#rasters for May and June.

generateMeanTemperatureRaster <- function(year){
  zPath <- paste("data/Temperature/",year,"/Raw/PRISM_tmean_stable_4kmM3_",year,"_all_bil.zip",sep="")
  if (file.exists(zPath) == F){
    stop(paste("Raw Temperature data not found. Please download from 'https://www.prism.oregonstate.edu/recent/' and place zip file in 'data/Temperature",year,"Raw'.",sep="/"))
  }
  months <- c("05","06")
  exts <- c("bil","bil.aux.xml","hdr","prj","stx")
  cores <- paste("PRISM_tmean_stable_4kmM3_",year,months,"_bil.",sep="")
  fileWhiteList <- paste(rep(cores,each=length(exts)),exts,sep="")
  xdir <- paste("data/Temperature",year,"Raw/temp",sep="/")
  unzip(zPath,fileWhiteList,exdir=xdir)
  may <- raster::raster(paste(xdir,"/",cores[1],"bil",sep=""))
  june <- raster::raster(paste(xdir,"/",cores[2],"bil",sep=""))
  mean <- raster::overlay(may,june,fun=function(a,b){return((a+b)/2)})
  unlink(xdir,recursive = T)
  sPath <- paste("data/Temperature",year,"mayJuneMean.grd",sep="/")
  recDir(sPath,wd=getwd())
  raster::writeRaster(mean,sPath)
}

#The sites are then assigned a mean temperature value by sampling the summer
#temperature raster generated above. The mean of the 100m surrounding the site
#is used. This means of attaching temperature data to BBS sites is identical to
#the approach used by Haddou et al. (2022)

generateTemperaturedSiteList <- function(year){
  print("Attaching temperatures to sites...")
  rPath <- paste("data/BBS",year,"sites_no_temperature.rda",sep="/")
  sites <- loadOrGenerate(rPath,generateSiteList,year=year)
  
  sites$Temperature <- rep(NA,length(sites$UID))
  tPath <- paste("data/Temperature",year,"mayJuneMean.grd",sep="/")
  tempRast <- loadOrGenerate(tPath,generateMeanTemperatureRaster,year=year,loadFunction=raster::raster)
  for (i in 1:length(sites$Temperature)){
    buf <- st_buffer(sites$Position[i],100)
    meanVal <- exact_extract(tempRast,buf,"mean")
    sites$Temperature[i] <- meanVal
  }
  sPath <- paste("data/BBS",year,"sites_complete.rda",sep="/")
  recDir(sPath,wd=getwd())
  save(sites, file=sPath)
}

#The final and by far the most dangerous function in this file deals with
#converting the BBS dataset into 8 binary rasters that Scalescape can interact
#with. It can easily take over a day to complete, and requires over 512GB of
#space to perform the process, and over 131GB to store the result. Additionally,
#it is the code in the project which I was least able to check. Proceed with
#caution!

generateBinaryRasters <- function(year){
  lPath <- paste(getwd(),"/data/",year,"/landcover/NLCD_",substr(year,3,4),"_LC.img")
  
  if (file.exists(lPath) == F){
    stop(paste("NLCD cover data not found. Please download from 'https://www.mrlc.gov/data?f%5B0%5D=category%3ALand%20Cover&f%5B1%5D=region%3Aconus' and place zip file in '",lPath,"'.",sep=""))
  }
  
  lcData <- raster::raster(lPath)
  attTab <- lcData@data@attributes[[1]]
  
  #Create a hash table, relating the number used for each cover type with its
  #displayed name.
  idLookup <- new.env()
  for (level in levels(attTab$NLCD.Land.Cover.Class)){
    if (level != ""){
      idLookup[[level]] <- attTab[attTab$NLCD.Land.Cover.Class == level,]$ID[[1]]
    }
  }
  
  #Create another has table, explaining which display names should match to our
  #combined cover types.
  lookup <- new.env()
  lookup[["Forest"]] <- c('Deciduous Forest','Evergreen Forest','Mixed Forest')
  lookup[["Grassland"]] <- c('Shrub/Scrub','Hay/Pasture','Herbaceous')
  lookup[["Cropland"]] <- c('Cultivated Crops')
  lookup[["Wetland"]] <- c('Woody Wetlands','Emergent Herbaceous Wetlands')
  lookup[["UrbanOpen"]] <- c('Developed, Open Space')
  lookup[["UrbanLow"]] <- c('Developed, Low Intensity')
  lookup[["UrbanMedium"]] <- c('Developed, Medium Intensity')
  lookup[["UrbanHigh"]] <- c('Developed, High Intensity')
  
  #These two hashes are used to create two lists, which are used by raster to
  #reclassify the cover data down to our eight groupings.
  
  ids <- c()
  tos <- c()
  i <- 1
  for (to in ls(lookup)){
    for (from in lookup[[to]]){
      ids <- append(ids,idLookup[[from]])
      tos <- append(tos,i)
    }
    i <- i +1
  }
  repF <- data.frame("ID"=ids,"To"=tos)
  
  reduced <- subs(lcData, repF)
  
  outPath <- paste("data",year,"bitmaps.gri",sep="/")
  #this is also untested!
  recDir(outPath,wd=getwd())
  pabBrick <- layerize(reduced,filename=outPath,datatype="LOG1s",options="COMPRESS=LZW")
}