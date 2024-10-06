#'
#'  SmoothHeat.R
#'
#'   Nadaraya-Watson style smooth regression using diffusion
#'
#'   Copyright (C) 2018-2024 Adrian Baddeley, Tilman Davies and Suman Rakshit
#' 
#'   $Revision: 1.3 $ $Date: 2024/10/06 01:26:29 $

SmoothHeat <- function(X, ...) {
  UseMethod("SmoothHeat")
}

SmoothHeat.im <- function(X, sigma, ...) {
  blurHeat(X, sigma, ...)
}

SmoothHeat.ppp <- function(X, sigma, ..., weights=NULL) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  marx <- marks(X)
  if(!is.vector(marx)) stop("Marks of X should be a numeric vector")
  marx <- as.numeric(marx)
  if(is.null(weights)) {
    numwt <- marx
    denwt <- NULL
  } else {
    check.nvector(weights, npoints(X), oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, npoints(X))
    numwt <- marx * weights
    denwt <- weights
  }
  Y <- unmark(X)
  numer <- densityHeat(Y, sigma, weights=numwt, ...)
  denom <- densityHeat(Y, sigma, weights=denwt, ...)
  return(numer/denom)
}

bw.SmoothHeatppp <- function(X, ..., weights=NULL,
                          srange=NULL, ns=16, sigma=NULL,
                          leaveoneout=TRUE, verbose=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  marx <- marks(X)
  if(!is.vector(marx)) stop("Marks of X should be a numeric vector")
  marx <- as.numeric(marx)
  if(is.null(weights)) {
    numwt <- marx
    denwt <- NULL
  } else {
    check.nvector(weights, npoints(X), oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, npoints(X))
    numwt <- marx * weights
    denwt <- weights
  }
  #' compute weighted and unweighted intensity estimates
  U <- unmark(X)
  aNumer <- HeatEstimates.ppp(U, ..., weights=numwt,
                              srange=srange, ns=ns, sigma=sigma,
                              leaveoneout=leaveoneout, verbose=verbose)
  aDenom <- HeatEstimates.ppp(U, ..., weights=denwt,
                              srange=srange, ns=ns, sigma=sigma,
                              leaveoneout=leaveoneout, verbose=verbose)
  h     <- aDenom$h
  hname <- aDenom$hname
  #' compute smoother
  zhat  <- aNumer$lambda/aDenom$lambda
  #' compute least squares cross-validation criterion
  zobs <- matrix(marx, nrow(zhat), ncol(zhat), byrow=TRUE)
  CV <- rowSums((zhat - zobs)^2)
  iopt <- which.min(CV)
  result <- bw.optim(CV, h, iopt,
                     criterion="Least squares cross-validation",
                     hname=hname,
                     unitname=unitname(X))
  return(result)
}

