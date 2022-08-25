#Returns a mapping of hand-picked colours for each cover type. Generally these
#are based on the colours used by the USGS, with some adjustments made for
#visibility/aesthetics.

getColEnv <- function(){
  colEnv <- new.env()
  colEnv[['Cropland']] <- rgb(202/255,145/255,70/255,1)
  colEnv[['Forest']] <- rgb(56/255,129/255,78/255,1)
  colEnv[['Grassland']] <- rgb(253/255,222/255,133/255,1)
  colEnv[['UrbanHigh']] <- rgb(181/255,0,0,1)
  colEnv[['UrbanMedium']] <- rgb(1,0,0,1)
  colEnv[['UrbanLow']] <- rgb(226/255,158/255,140/255,1)
  colEnv[['UrbanOpen']] <- rgb(232/255,209/255,209/255,1) 
  colEnv[['Wetland']] <- rgb(0,0.5,0.5,1)
  return(colEnv)
}

#Sets up par to optimally display an arbitrary number of square plots.

parTile <- function(elems){
  se <- ceiling(sqrt(elems))
  par(mfcol=c(se,ceiling(elems/se)),pty="s",mai=c(0.1,0.1,0.1,0.1))
}

#Naive load balancer. divides into n even chunks, then adds remainder to the
#last one. This leads to a maximum discrepency of bins-1, which is negligable, 
#since this is mostly to help with load balancing.

naiveBalance <- function(num,bins){
  return(append(rep(num%/%bins,bins-1),(num%/%bins)+num%%bins))
}

#For a given vector, returns a list of vectors identical except for one item
#being removed.

dropOnePerms <- function(vec){
  perms <- list()
  for (i in 1:length(vec)){
    perms[[i]] <- append(head(vec,i-1),tail(vec,length(vec)-i))
  }
  return(perms)
}

#Searches a vectoroid for a given value, return the index of the first find or
#-1 iff absent.

findIn <- function(val,vec){
  if (length(vec) == 0){
    return(-1)
  }
  for (i in 1:length(vec)){
    if (vec[i] == val){
      return(i)
    }
  }
  return(-1)
}

#Takes a path and combines it with the given working directory, then creates
#directories as necessary until the path is complete. Default behaviour
#(overidden by the noFile flag) is to ignore the last segment of the path. This
#means that handing the path you wish to save a file to will automatically
#create the directory nest it should be in.

recDir <- function(path,wd="",noFile=F){
  dirTree <- strsplit(path,"/")[[1]]
  if (noFile == F){
    dirTree <- head(dirTree,-1) 
  }
  if (wd == ""){
    wd <- getwd()
  }
  acc <- wd
  for (i in 1:length(dirTree)){
    acc <- paste(acc,dirTree[i],sep="/")
    if (dir.exists(acc) == F){
      dir.create(acc)
    }
  }
}

#This is the core of the lazy-evaluating "load-or-generate" system used by this
#project's code. loadOrGenerate takes a filepath "path" and a generator function
#"genfun". If the file exists, it simply loads and returns it. If not, it will
#call the generator function, passing along any additional arguments. Generator
#functions will always save their output to disk, so once the function has run,
#loadOrGenerate will load from disk as usual. By default, "load" will be used to
#read the saved file, but another function can be handed in instead using the
#"loadFunction" parameter.

loadOrGenerate <- function(path,genFun,loadFunction=NULL,...){
  if (file.exists(path) == F){
    recDir(path)
    genFun(...)
  }
  if (is.null(loadFunction)){
    loadFunction <- function(ppath){
      obj <- get(load(ppath)[1])
      return(obj)
    }
  }
  return(loadFunction(path))
}