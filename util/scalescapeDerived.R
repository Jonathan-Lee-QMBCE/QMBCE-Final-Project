# All the functions in this file were either directly copied from the
# 'Scalescape' package (https://github.com/benjaminiuliano/scalescape), or are
# based so heavily on its code as to be derivative. The package is licensed with
# with a GPL 3 license, and is therefore permissible to edit, but these
# functions were quarantined from non-derivative code for clarity.
#
# This code is not very well commented, since I honestly do not understand the
# subtleties of it, and the original developers broadly didn't comment their
# code.

#Modified version of the plot function that actually works. The unmodified
#version uses x in two contexts at once during curve construction, causing a
#fatal crash.
modded.plot.scalescape <- function(x) {
  if(!is.element("scalescape.boot", class(x))) {
    sol <- modded.weighting_funct(x$opt.range, mod0=x$mod0, landscape.formula = x$landscape.formula,
                                  data = x$data, max.Dist = x$max.Dist, landscape.vars=x$landscape.vars, weight.fn = x$weight.fn, return.coef = T)
    par(mfrow=c(length(x$opt.range),2), mai=c(0.1,0.1,.15,.15))
    for(i.range in 1:length(x$opt.range)){
      w <- data.frame(range=x$max.Dist[i.range] * .005*(1:199), logLik = NA, b = NA)
      for(i in 1:199) {
        opt.par <- x$opt.range
        opt.par[i.range] <- .005*i
        z <- modded.weighting_funct(par=opt.par, mod0=x$mod0, landscape.formula = x$landscape.formula,
                                    data = x$data, max.Dist = x$max.Dist, landscape.vars=x$landscape.vars, weight.fn = x$weight.fn, return.coef = T)
        if(length(z)>1){
          w$logLik[i] <- z$logLik
          #WHY 3?!
          w$b[i] <- z$coef[3]
        }
      }
      #1
      plot(logLik ~ range, data=w, typ="l", xlab="Distance", ylab="logLik")
      points(x$opt.range[i.range], x$logLik[i.range], col="red")
      mtext(side=3, paste0("Range for ",names(x$opt.range[i.range])," = ",round(x$opt.range[i.range], digits=3)))
      
      #2
      # Gaussian weightings
      
      if(x$weight.fn=="Gaussian"){ 
        #Commented out below is the original broken version.
        #curve(exp(-0.5*(x/(x$opt.range[i.range]))^2), from = 0, to = x$max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
        optr <- x$opt.range[i.range]
        curve(exp(-0.5*(x/(optr))^2), from = 0, to = x$max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
        abline(v=x$opt.range[i.range], col="red", lty=2)
      }
      # exponential weightings
      #not tested or fixed, but I imagine the same problems exist
      if(x$weight.fn=="exponential") curve(exp(-x/(x$opt.range[i.range])),
                                           from = 0, to = x$max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=x$opt.range[i.range], col="red", lty=2)
    }
  }else{
    par(mfcol=c(1,1 + length(x$mod.full$opt.range)), mai=c(0.3,0.3,.15,.15))
    
    #1
    if(length(x$mod.full$opt.range) == 1){
      hist(x$mod.full$max.Dist * x$logLik.values$opt.range.boot,
           xlab="Range", main=paste0("Bootstrap for ", names(x$mod.full$opt.range)))
    }else{
      for(i.range in 1:length(x$mod.full$opt.range))
        hist(x$mod.full$max.Dist[i.range] * x$logLik.values[,paste0("opt.range.boot.",i.range)],
             xlab="Range", main=paste0("Bootstrap for ", names(x$mod.full$opt.range)[i.range]))
    }
    
    #2
    hist(x$logLik.values$dev, xlab="Deviance", main=paste0("Bootstrap P = ", round(x$P, digits=4)),
         xlim=c(-.01,max(c(x$logLik.values$dev,x$logLik.values$dev))))
    lines(c(x$logLik.values$dev,x$logLik.values$dev), c(0,nrow(x$logLik.values)), col="red")
  }
}

#Reduced plotting with just the curves. The likelihood plots are mostly
#nonsense anyway, and take a very long time to generate.
plotScalescapeCurves <- function(x, ...) {
  if(!is.element("scalescape.boot", class(x))) {
    if (length(x$opt.range) == 1){
      par(mfrow=c(1,1), mai=c(0.5,0.5,0.5,0.5))  
    }
    else{
      par(mfrow=c((length(x$opt.range)%/%2)+(length(x$opt.range)%%2),2), mai=c(0.3,0.3,.3,.3))  
    }
    for(i.range in 1:length(x$opt.range)){
      if(x$weight.fn=="Gaussian"){ 
        optr <- x$opt.range[i.range]
        curve(exp(-0.5*(x/(optr))^2), from = 0, to = x$max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
        abline(v=x$opt.range[i.range], col="red", lty=2)
      }
      # exponential weightings
      if(x$weight.fn=="exponential") curve(exp(-x/(x$opt.range[i.range])),
                                           from = 0, to = x$max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=x$opt.range[i.range], col="red", lty=2)
      mtext(side=3, names(x$opt.range)[i.range])
    }
  }
}

assessLandscapeData <- function(landDat,weight.fn,parVal,max.Dist){
  if(is.element("matrix",is(landDat))) {
    Dist <- landDat[,1]
    if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/parVal/max.Dist)^2)
    if(weight.fn=="exponential") weighting <- exp(-Dist/parVal/max.Dist)
    
    return(t(t(weighting) %*% landDat[, -1]/sum(weighting)))
  }
  if(is.element("list",is(landDat))) {
    landscape.list <- landDat
    Dist <- landscape.list[[1]][,1]
    if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/parVal/max.Dist)^2)
    if(weight.fn=="exponential") weighting <- exp(-Dist/parVal/max.Dist)
    weighting <- as.matrix(weighting)/mean(weighting)
    
    X.landclass.sub <- NULL
    for(i.list in 1:length(landscape.list)){
      X.landclass.sub <- cbind(X.landclass.sub, as.matrix(t(t(weighting) %*% landscape.list[[i.list]][, -1]/sum(weighting)), ncol=1))
    }
    if(!is.null(names(landscape.list))) colnames(X.landclass.sub) <- names(landscape.list)
    return(X.landclass.sub)
    
  }
}

#Modified so that the landscape variables can be calculated in parallel,
#although if I recall correctly, this actually made things worse. Regardless,
#some version of the weighting function is needed here, since the one used
#within the Scalescape package is not accessible. 
modded.weighting_funct <- function(par, mod0, landscape.formula,
                                   data, landscape.vars, max.Dist, weight.fn = "Gaussian", return.coef = F, n.cores=1){
  
  new.vars <- names(landscape.vars)
  mod0.vars <- row.names(attr(terms(formula(mod0)), "factors"))
  mod.landscape.vars <- c(row.names(attr(terms(formula(mod0)), "factors")), new.vars)
  
  
  landscape.data <- data
  i.count <- 0
  
  if (n.cores < 0){
    lDList <- list()
    for (i in 1:length(new.vars)){
      lDList[[i]] <- landscape.vars[[new.vars[i]]]
    }
    cl <- makeCluster(n.cores)
    landScapeData <- clusterMap(cl=cl,max.Dist=max.Dist,parVal=par,landDat=lDList,fun=assessLandscapeData,MoreArgs=list(weight.fn=weight.fn))
    stopCluster(cl) 
    for (i in 1:length(new.vars)){
      landscape.data[as.character(new.vars[i])] <- landScapeData[[i]]
    }
  }
  else{
    for(i.terms in new.vars){
      if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
        i.count <- i.count + 1
        
        Dist <- landscape.vars[[i.terms]][,1]
        if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
        if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])
        
        landscape.data[as.character(i.terms)] <- t(t(weighting) %*% landscape.vars[[i.terms]][, -1]/sum(weighting))
      }
      if(is.element("list",is(landscape.vars[[i.terms]]))) {
        i.count <- i.count + 1
        landscape.list <- landscape.vars[[i.terms]]
        
        Dist <- landscape.list[[1]][,1]
        if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
        if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])
        weighting <- as.matrix(weighting)/mean(weighting)
        
        X.landclass.sub <- NULL
        for(i.list in 1:length(landscape.list)){
          X.landclass.sub <- cbind(X.landclass.sub, as.matrix(t(t(weighting) %*% landscape.list[[i.list]][, -1]/sum(weighting)), ncol=1))
        }
        if(!is.null(names(landscape.list))) colnames(X.landclass.sub) <- names(landscape.list)
        landscape.data[as.character(i.terms)] <- X.landclass.sub
        
      }
    }
  }
  formula <- formula(mod0)
  # show(par)
  mod <- try(update(mod0, landscape.formula, data=landscape.data))
  
  # if(is(mod)[1] == "try-error") return(-10^10)
  if(return.coef){
    if(!is.element(class(mod0)[1], c("lmerMod", "glmerMod"))) coef <- mod$coef else coef <- summary(mod)$coef[,1]
    return(list(logLik=logLik(mod), coef = coef, mod = mod, data.with.landscape = landscape.data))
  }else{
    return(-logLik(mod))
  }
}

#Uses a fitted scalescape model to assign values to landscape variables for new
#data. This functionality is not possible with the Scalescape library, a rather
#glaring ommission, since it prevents you from using a fitted Scalescape model
#to make predictions.
getWeighted <- function(par, mod0, landscape.formula, data, landscape.vars, max.Dist, weight.fn = "Gaussian"){
  new.vars <- names(landscape.vars)
  mod0.vars <- row.names(attr(terms(formula(mod0)), "factors"))
  mod.landscape.vars <- c(row.names(attr(terms(formula(mod0)), "factors")), new.vars)
  
  landscape.data <- data
  i.count <- 0
  for(i.terms in new.vars){
    if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1
      
      Dist <- landscape.vars[[i.terms]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])
      
      landscape.data[as.character(i.terms)] <- t(t(weighting) %*% landscape.vars[[i.terms]][, -1]/sum(weighting))
    }
    if(is.element("list",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1
      landscape.list <- landscape.vars[[i.terms]]
      
      Dist <- landscape.list[[1]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])
      weighting <- as.matrix(weighting)/mean(weighting)
      
      X.landclass.sub <- NULL
      for(i.list in 1:length(landscape.list)){
        X.landclass.sub <- cbind(X.landclass.sub, as.matrix(t(t(weighting) %*% landscape.list[[i.list]][, -1]/sum(weighting)), ncol=1))
      }
      if(!is.null(names(landscape.list))) colnames(X.landclass.sub) <- names(landscape.list)
      landscape.data[as.character(i.terms)] <- X.landclass.sub
      
    }
  }
  return(landscape.data)
}

#Not actually modified at all, but a local copy is needed since the one within
#the Scalescape package is inaccessible.
modded.dist_weight <- function(mod0, landscape.vars, landscape.formula,
                        data = NULL, weight.fn = "Gaussian",  plot.fits = TRUE,
                        lower = NULL, upper = NULL, init.range = NULL, n.partition = 10,
                        opt.range = NULL, optim.method = "L-BFGS-B", n.cores=1){
  new.vars <- names(landscape.vars)
  mod0.vars <- row.names(attr(terms(formula(mod0)), "factors"))
  mod.landscape.vars <- c(row.names(attr(terms(formula(mod0)), "factors")), new.vars)
  
  if(!is.null(opt.range) & (length(opt.range) != length(new.vars))) stop("If you specify opt.range, the number of values must equal the number of new spatial elements estimated.")
  if(!is.null(init.range) & (length(init.range) != length(new.vars))) stop("If you specify init.range, the number of values must equal the number of new spatial elements estimated.")
  
  if(is.null(data) & is.element(class(mod0)[1], c("lm","glm"))) data <- mod0$model
  if(is.null(data) & is.element(class(mod0)[1], c("lmerMod", "glmerMod"))) data <- model.frame(mod0)
  if(is.null(data) & class(mod0)[1] == "lme") data <- mod0$data
  if(is.null(data) & class(mod0)[1] == "gls"){
    if(is.null(data)) stop("For gls() you need to specify data")
    data <- data
  }
  if(class(mod0)[1] == "gls" && mod0$method == "REML"){
    warning("For gls(), fitting show be with method = ML. Therefore, we refit your model.")
    mod0 <- stats::update(mod0, method="ML")
  }
  
  max.Dist <- NULL
  for(i.terms in new.vars){
    if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
      Dist <- landscape.vars[[i.terms]][,1]
      max.Dist <- c(max.Dist, max(Dist))
    }
    if(is.element("list",is(landscape.vars[[i.terms]]))) {
      landscape.list <- landscape.vars[[i.terms]]
      Dist <- landscape.list[[1]][,1]
      max.Dist <- c(max.Dist, max(Dist))
    }
  }
  
  if(is.null(opt.range)){
    if(is.null(init.range)) {
      init.range <- rep(0.2, length(new.vars))
      opt.init.range <- init.range
      z <- modded.weighting_funct(par=init.range, mod0=mod0, landscape.formula = landscape.formula, data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T, n.cores=n.cores)
      opt.logLik <- z$logLik
      for(i.range in 1:length(new.vars)){
        for(i in 1:n.partition) {
          par <- opt.init.range
          par[i.range] <- i/(n.partition+1)
          z <- modded.weighting_funct(par=par, mod0=mod0, landscape.formula = landscape.formula,
                               data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T, n.cores)
          
          if(z$logLik > opt.logLik){
            opt.logLik <- z$logLik
            opt.init.range <- par
          }
          # show(c(z$logLik, par))
        }
      }
      init.range <- opt.init.range
    }else{
      init.range <- init.range/max.Dist
    }
    # if(length(init.range)==1){
    # optim.method="Brent"
    # }else{
    #	optim.method="L-BFGS-B"
    # }
    
    if(is.element(optim.method, c("Brent", "L-BFGS-B"))){
      if(is.null(lower)) {
        lower=rep(.01, length(new.vars))
      }else{
        lower=lower/max.Dist
      }
      if(is.null(upper)) {
        upper=rep(.99, length(new.vars))
      }else{
        upper=upper/max.Dist
      }
      
    }else{
      lower=NULL
      upper=NULL
    }
    if (n.cores == 1 || optim.method != "L-BFGS-B"){
      print(paste("using",n.cores,"cores."))
      opt.range <- optim(par=init.range, fn=modded.weighting_funct, method=optim.method, upper = upper, lower = lower, mod0=mod0, data = data, max.Dist = max.Dist, landscape.formula = landscape.formula, landscape.vars=landscape.vars, weight.fn = weight.fn, n.cores=n.cores)$par
    }
    else{
      print("optimising with multicore method.")
      cl <- makeCluster(n.cores)
      clusterEvalQ(cl, library("scalescape"))
      clusterEvalQ(cl, library("lme4"))
      opt.range <- optimParallel(par=init.range, fn=modded.weighting_funct, parallel=list(cl=cl), method=optim.method, upper = upper, lower = lower, mod0=mod0, data = data, max.Dist = max.Dist, landscape.formula = landscape.formula, landscape.vars=landscape.vars, weight.fn = weight.fn)$par
      stopCluster(cl) 
    }
    names(opt.range) <- new.vars
  } else {
    opt.range <- opt.range/max.Dist
  }
  sol <- modded.weighting_funct(opt.range, mod0 = mod0, landscape.formula = landscape.formula,
                         data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T,n.cores=n.cores)
  
  npar <- length(sol$coef) + length(opt.range) + 1
  AIC <- 2*npar - 2*sol$logLik
  BIC <- 2*npar*log(nrow(data)) - 2*sol$logLik
  
  if(plot.fits){
    par(mfrow=c(length(opt.range),2), mai=c(1,1,.3,.3))
    for(i.range in 1:length(opt.range)){
      w <- data.frame(range=max.Dist[i.range] * .005*(1:199), logLik = NA, b = NA)
      for(i in 1:199) {
        opt.par <- opt.range
        opt.par[i.range] <- .005*i
        z <- modded.weighting_funct(par=opt.par, mod0=mod0, landscape.formula = landscape.formula,
                             data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T,n.cores=n.cores)
        if(length(z)>1){
          w$logLik[i] <- z$logLik
          w$b[i] <- z$coef[3]
        }
      }
      
      #1
      plot(logLik ~ range, data=w, typ="l", xlab="Distance", ylab="logLik")
      points(max.Dist[i.range] * opt.range[i.range], sol$logLik, col="red")
      mtext(side=3, paste0("Range for ",new.vars[i.range]," = ",round(max.Dist[i.range]*opt.range[i.range], digits=3)))
      
      #2
      # Gaussian weightings
      if(weight.fn=="Gaussian") curve(exp(-0.5*(x/(max.Dist[i.range] * opt.range[i.range]))^2),
                                      from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)
      
      # exponential weightings
      if(weight.fn=="exponential") curve(exp(-x/(max.Dist[i.range] * opt.range[i.range])),
                                         from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)
    }
  }
  
  to.return <- list(opt.range = max.Dist * opt.range,
                    logLik = sol$logLik,
                    AIC = AIC,
                    BIC = BIC,
                    npar = npar,
                    data = data,
                    coef = sol$coef,
                    mod = sol$mod,
                    mod0 = mod0,
                    landscape.formula = landscape.formula,
                    data.with.landscape = sol$data.with.landscape,
                    landscape.vars = landscape.vars,
                    weight.fn = weight.fn,
                    max.Dist = max.Dist)
  class(to.return) <- "scalescape"
  return(to.return)
}