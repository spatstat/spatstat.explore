#'
#'   relriskHeat.R
#'
#'   Relative risk/conditional probability using diffusion smoothing
#'
#'   Copyright (C) 2018-2024 Adrian Baddeley, Tilman Davies and Suman Rakshit
#'
#'   $Revision: 1.7 $ $Date: 2024/10/06 02:55:11 $
#'
#' 

relriskHeat <- function(X,...) {
  UseMethod("relriskHeat")
}

relriskHeat.ppp <- function(X,..., sigmaX=NULL, weights=NULL){
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  nX <- npoints(X)
  
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  
  type <- marks(X)
  
  if(length(sigmaX)) {
    check.nvector(sigmaX, nX)
    sigmaX <- split(sigmaX, type)
  } else sigmaX <- rep(list(NULL), ntypes)

  if(length(weights)) {
    check.nvector(weights, nX)
    weights <- split(weights, type)
  } else weights <- rep(list(NULL), ntypes)

  Deach <- mapply(densityHeat, x=Y, sigmaX=sigmaX, weights=weights,
                  MoreArgs=list(...), SIMPLIFY=FALSE)
  
  Dall <- Reduce("+", Deach)
  
  probs <- solapply(Deach, "/", e2=Dall)
  
  return(probs)
}

bw.relriskHeatppp <- function(X, ...,
                              method=c("likelihood", "leastsquares"),
                              weights=NULL,
                              srange=NULL, ns=16, sigma=NULL,
                              leaveoneout=TRUE, verbose=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  method <- match.arg(method)
  U <- unmark(X)
  Y <- split(X)
  Denominator <- HeatEstimates.ppp(U, ..., weights=weights,
                                   srange=srange, ns=ns, sigma=sigma,
                                   leaveoneout=leaveoneout, verbose=verbose)
  h     <- Denominator$h
  hname <- Denominator$hname
  # extract denominator value for each sigma (row) and each data point (col)
  lambda.denom <- Denominator$lambda
  #'
  if(is.null(weights)) {
    Numerators <- lapply(Y, HeatEstimates.ppp, ...,
                         srange=srange, ns=ns, sigma=sigma,
                         leaveoneout=leaveoneout, verbose=verbose)
  } else {
    check.nvector(weights, npoints(X), oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, npoints(X))
    wsplit <- split(weights, marks(X))
    Numerators <- mapply(HeatEstimates.ppp,
                         X=Y, weights=wsplit,
                         MoreArgs = list(...,
                                         srange=srange, ns=ns, sigma=sigma,
                                         leaveoneout=leaveoneout,
                                         verbose=verbose),
                         SIMPLIFY=FALSE)
  }
  #' extract estimates of numerator terms
  lamlist <- lapply(Numerators, getElement, name="lambda")
  #' reassemble into original position
  #' (tried to do this with 'unsplit' but it's too messy)
  opos <- split(seq_len(npoints(X)), marks(X))
  lambda.numer <- matrix(, nrow=nrow(lambda.denom), ncol=ncol(lambda.denom))
  for(k in seq_along(opos)) {
    if(length(opos.k <- opos[[k]]))
      lambda.numer[ , opos.k] <- lamlist[[k]]
  }
  #' compute predicted probability of observations
  phat  <- lambda.numer/lambda.denom
  #' compute cross-validation criterion
  switch(method,
         likelihood = {
           CV <- -rowMeans(log(phat))
           cname <- "Likelihood cross-validation"
         },
         leastsquares = {
           CV <- rowMeans((1 - phat)^2)
           cname <- "Least squares cross-validation"
         })
  result <- bw.optim(CV, h, criterion=cname, hname=hname,
                     unitname=unitname(X))
  return(result)
}

