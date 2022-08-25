library(dplyr)

#works on assumption that if you move stop every 50th point along
polyLineToPoints <- function(entry){
  stepLength <- (entry$RTELENG/48)/1000
  pointMat <- as.matrix(routeShapes$geometry[[1]][])
  pointMat[ , c(1,2)]  <- pointMat[ , c(2,1)]
  stops <- data.frame("lat"=pointMat[1,2],"lon"=pointMat[1,1],"dis"=0)
  acc <- 0
  for (i in 2:length(pointMat[,1])){
    print(i)
    ndis <- spDistsN1(pointMat[(i-1):i,],pointMat[i-1,],T)[2]
    if (ndis+acc >= stepLength || i == length(pointMat[,1])){
      print("overspill")
      overspill <- (ndis+acc-stepLength)
      print(ndis+acc)
      acc = overspill
      pointLat <- pointMat[i-1,2] + (pointMat[i,2]-pointMat[i-1,2])*(1-overspill)
      pointLon <- pointMat[i-1,1] + (pointMat[i,1]-pointMat[i-1,1])*(1-overspill)
      #pointLat <- pointMat[i,2]
      #pointLon <- pointMat[i,1]
      
      tMat <- matrix(c(stops$lon[length(stops$lon)],stops$lat[length(stops$lon)],pointLon,pointLat),ncol=2)
      print(tMat)
      dfl <- spDistsN1(tMat,tMat[1,],T)[2]
      stops <- rbind(stops,data.frame("lat"=pointLat,"lon"=pointLon,"dis"=dfl))
    }
    else{
      acc <- acc+ndis
    }
  }
  return(stops)
}