#' monotonicRhoFun.R
#'
#' monotonic regression of (intensity of) point pattern X on covariate Z
#'
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner


monotonicRhoFun <- function(X, Z, increasing=FALSE,
                      ..., weights=NULL, subset=NULL, baseline=NULL) {
  stopifnot(is.ppp(X) || is.lpp(X))
  trap.extra.arguments(..., .Context="In function monotonicRhoFun()")
  if(is.null(weights)) {
    weights <- 1
  } else {
    check.nvector(weights, npoints(X), oneok=TRUE)
    stopifnot(all(weights >= 0))
  }
  if(!is.ppp(baseline)) {
    #' usual case
    nullmodel <- resolveNullModel(baseline, X)
    stuff <- spatialCovariateEvidence(nullmodel, Z,
                                      subset=subset, jitter=FALSE)$values
    lambda <- monotonicRhoFunCalc(x = stuff$ZX,
                            massx = weights,
                            z = stuff$Zvalues,
                            weightz = stuff$weights,
                            increasing = increasing)
  } else {
    #' baseline is a point pattern
    Y <- baseline
    Z <- digestCovariates(Z, W=Window(X))[[1]]
    ZX <- evaluateCovariate(Z, X)
    ZY <- evaluateCovariate(Z, Y)
    lambda <- monotonicRhoFunCalc(x       = ZX,
                            massx   = weights,
                            z       = ZY,
                            weightz = 1,
                            increasing = increasing)
  }
  attr(lambda, "call") <- sys.call()
  return(lambda)
}

monotonicRhoFunCalc <- function(x, z, massx=1, weightz=1, increasing=FALSE) {
  #' x = observed values
  #' massx = optional weights for x  (usually 1)
  #' z = sample from reference distribution
  #' weightz = weights for z (eg pixel area)
  if(length(massx) == 0) massx <- 1
  if(length(massx) == 1) massx <- rep(massx, length(x))
  if(length(weightz) == 0) weightz <- 1
  if(length(weightz) == 1) weightz <- rep(weightz, length(z))
  #' sort x into increasing order
  ox <- order(x)
  x <- x[ox]
  massx <- massx[ox]
  #' cdf of reference distribution
  g <- ewcdf(z, weightz, normalise=FALSE)
  #'
  areas <- g(x)
  if(increasing) areas <- sum(weightz) - rev(areas)
  ## maximum upper sets algorithm
  y <- numeric(0)
  a <- massx
  b <- diff(c(0, areas))
  while(length(b) > 0) {
    u <- cumsum(a)/cumsum(b)
    if(any(bad <- !is.finite(u))) # divide by zero etc
      u[bad] <- max(u[!bad], 0)
    k <- which.max(u)
    y <- c(y, rep(u[k], k))
    a <- a[-(1:k)]
    b <- b[-(1:k)]
  }
  y <- c(y, 0)
  if(increasing) y <- rev(y)
  lambda <- stepfun(x = x, y=y, right=TRUE, f=1)
  return(lambda)
}

