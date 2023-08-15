##
## auc.R
##
##  Calculate ROC curve or area under it
##
## $Revision: 1.17 $ $Date: 2023/08/15 07:44:11 $

roc <- function(X, ...) { UseMethod("roc") }

roc.ppp <- function(X, covariate, ..., high=TRUE) {
  nullmodel <- exactppm(X)
  result <- rocData(covariate, nullmodel, ..., high=high)
  return(result)
}

rocData <- function(covariate, nullmodel, ...,
                    high=TRUE,
                    p=seq(0, 1, length=1024)) {
  d <- spatialCDFframe(nullmodel, covariate, ...)
  U <- d$values$U
  ec <- if(high) ecdf(1-U) else ecdf(U)
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
  }
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

rocModel <- function(lambda, nullmodel, ..., high,
                     p=seq(0, 1, length=1024)) {
  if(!missing(high))
    warning("Argument 'high' is ignored when computing ROC for a fitted model")
  d<- spatialCDFframe(nullmodel, lambda, ...) 
  U <- d$values$U
  ec <- ecdf(1-U) 
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
  }
  fobs <- ec(p)
  FZ <- d$values$FZ
  FZinverse <- quantilefun.ewcdf(FZ)
  lambdavalues <- if(is.im(lambda)) lambda[] else unlist(lapply(lambda, "["))
  F1Z <- ewcdf(lambdavalues, lambdavalues/sum(lambdavalues))    
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


