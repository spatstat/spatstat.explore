##
##
##  auc.R
##
##  Calculate Area Under ROC curve
##       (in spatstat.explore)
##
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner
##

auc <- function(X, ...) { UseMethod("auc") }

needROC <- function(...) {
  ## these arguments require the use of roc()
  any(c("observations", "baseline", "weights") %in% names(list(...)))
}

auc.ppp <- function(X, covariate, ...,
                    high=TRUE,  subset=NULL) {
  verifyclass(X, "ppp")
  if(needROC(...)) {
    ro <- roc(X, covariate, ..., high=high, subset=subset)
    result <- auc(ro)
  } else {
    nullmodel <- exactppm(X)
    result <- aucData(covariate, nullmodel, ..., high=high, subset=subset)
  }
  return(result)
}


aucData <- function(covariate, nullmodel, ..., high=TRUE,
                    interpolate=FALSE, jitter=FALSE, subset=NULL,
                    covariateAtPoints=discrimAtPoints,
                    discrimAtPoints=NULL) {
  d <- spatialCDFframe(nullmodel, covariate, ...,
                       covariateAtPoints=covariateAtPoints,
                       interpolate=interpolate, jitter=jitter,
                       subset=subset)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU)
  return(result)
}




auc.im <- function(X, covariate, ..., high=TRUE) {
  ro <- roc(X, covariate, ..., high=high)
  auc(ro)
}

auc.roc <- function(X, ...) {
  with(X, colMeans(.))
}


## AUC methods for other classes
auc.spatialCDFframe <- function(X, ..., high=TRUE) {
  trap.extra.arguments(...)
  U <- X$values$U
  EU <- mean(U)
  result <- if (high) EU else (1 - EU)
  return(result)
}

auc.bermantest <- function(X, ..., high=TRUE) {
  auc(X[["fram"]], high=high, ...)
}

auc.cdftest <- function(X, ..., high=TRUE) {
  auc(attr(X, "frame"), high=high, ...)
}
