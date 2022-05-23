##
## auc.R
##
##  Calculate ROC curve or area under it
##
## $Revision: 1.15 $ $Date: 2022/05/22 04:08:02 $

roc <- function(X, ...) { UseMethod("roc") }

roc.ppp <- function(X, covariate, ..., high=TRUE) {
  nullmodel <- exactppm(X)
  result <- rocData(covariate, nullmodel, ..., high=high)
  return(result)
}

rocData <- function(covariate, nullmodel, ..., high=TRUE) {
  d <- spatialCDFframe(nullmodel, covariate, ...)
  U <- d$values$U
  ec <- if(high) ecdf(1-U) else ecdf(U)
  p <- seq(0,1,length=1024)
  df <- data.frame(p=p, fobs=ec(p), fnull=p)
  result <- fv(df,
               argu="p",
               ylab=quote(roc(p)),
               valu="fobs",
               fmla= . ~ p,
               desc=c("fraction of area",
                      "observed fraction of points",
                      "expected fraction if no effect"),
               fname="roc")
  fvnames(result, ".") <- c("fobs", "fnull")
  return(result)
}

rocModel <- function(lambda, nullmodel, ..., high) {
  if(!missing(high))
    warning("Argument 'high' is ignored when computing ROC for a fitted model")
  d<- spatialCDFframe(nullmodel, lambda, ...) 
  U <- d$values$U
  ec <- ecdf(1-U) 
  p <- seq(0,1,length=1024)
  fobs <- ec(p)
  FZ <- d$values$FZ
  lambdavalues <- if(is.im(lambda)) lambda[] else unlist(lapply(lambda, "["))
  F1Z <- ewcdf(lambdavalues, lambdavalues/sum(lambdavalues))    
  pZ <- get("y", environment(FZ))
  qZ <- get("x", environment(FZ))
  FZinverse <- approxfun(pZ, qZ, rule=2)
  ftheo <- 1 - F1Z(FZinverse(1-p))
  df <- data.frame(p=p, fobs=fobs, ftheo=ftheo, fnull=p)
  result <- fv(df,
               argu="p",
               ylab=quote(roc(p)),
               valu="fobs",
               fmla = . ~ p,
               desc=c("fraction of area",
                 "observed fraction of points",
                 "expected fraction of points",
                 "expected fraction if no effect"),
               fname="roc")
  fvnames(result, ".") <- c("fobs", "ftheo", "fnull")
  return(result)
}


## Code for roc.ppm, roc.slrm, roc.kppm is moved to spatstat.model

#    ......................................................

auc <- function(X, ...) { UseMethod("auc") }

auc.ppp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(exactppm(X), covariate, ...)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU) 
  return(result)
}

## Code for auc.ppm, auc.slrm, auc.kppm is moved to spatstat.model


