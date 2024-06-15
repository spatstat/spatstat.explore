#'
#' Spatially weighted quantile
#'

SpatialMedian <- function(X, ...) {
  UseMethod("SpatialMedian")
}

SpatialQuantile <- function(X, prob=0.5, ...) {
  UseMethod("SpatialQuantile")
}

#' methods for 'ppp' class

SpatialMedian.ppp <- function(X, sigma=NULL, ...,
                             type=4,
                             at=c("pixels", "points"), leaveoneout=TRUE,
                             weights=NULL, 
                             edge=TRUE, diggle=FALSE, verbose=FALSE) {
  SpatialQuantile.ppp(X, sigma=sigma, prob=0.5,
                     ...,
                     type=type,
                     at=at, leaveoneout=leaveoneout,
                     weights=weights,
                     edge=edge, diggle=diggle, verbose=verbose)
}

SpatialQuantile.ppp <- function(X, prob=0.5, sigma=NULL, ...,
                               type=1,
                               at=c("pixels", "points"), leaveoneout=TRUE,
                               weights=NULL,
                               edge=TRUE, diggle=FALSE, verbose=FALSE) {
  if(!is.ppp(X)) stop("X should be a point pattern")
  if(!is.marked(X)) stop("The point pattern X should have marks")
  check.1.real(prob)
  stopifnot(prob >= 0)
  stopifnot(prob <= 1)
  at <- match.arg(at)
  atName <- switch(at, pixels="pixels", points="data points")
  check.1.integer(type)
  type <- as.integer(type)
  if(!any(type == c(1L,4L))) 
    stop(paste("Quantiles of type", type, "are not supported"), call.=FALSE)
  ## extract marks
  X <- coerce.marks.numeric(X)
  m <- marks(X)
  ## multiple columns of marks?
  if(!is.null(dim(m)) && ncol(m) > 1) {
    ## compute separately for each column
    Xlist <- unstack(X)
    Zlist <- lapply(Xlist, SpatialQuantile, prob=prob, sigma=sigma, ...,
                    type=type,
                    at=at, leaveoneout=leaveoneout, weights=weights,
                    edge=edge, diggle=diggle, verbose=verbose)
    ZZ <- switch(at,
                 pixels = as.imlist(Zlist),
                 points = do.call(data.frame, Zlist))
    return(ZZ)
  }
  ## single column of marks
  m <- as.numeric(m)
  nX <- npoints(X)
  #' unique mark values
  um <- sort(unique(m))
  Num <- length(um)
  #' trivial cases
  if(nX == 0 || ((Num == 1) && leaveoneout)) {
    Z <- switch(at,
                pixels = as.im(NA_real_, W=Window(X), ...),
                points = rep(NA_real_, nX))
    attr(Z, "sigma") <- sigma
    return(Z)
  }
  if(Num == 1) {
    Z <- switch(at,
                pixels = as.im(um[1], W=Window(X), ...),
                points = rep(um[1], nX))
    attr(Z, "sigma") <- sigma
    return(Z)
  }
  #' numerical weights 
  if(!is.null(weights)) {
    check.nvector(weights, nX, vname="weights")
    if(any(weights < 0)) stop("Negative weights are not permitted")
    if(sum(weights) < .Machine$double.eps)
      stop("Weights are numerically zero; quantiles are undefined", call.=FALSE)
  }
  #' start main calculation
  ## bandwidth selector
  if(is.function(sigma))
    sigma <- sigma(X, ...)
  #' edge correction has no effect if diggle=FALSE
  #'     (because uniform edge correction cancels)
  edge <- edge && diggle
  #' smoothed intensity of entire pattern
  UX <- unmark(X)
  LX <- density(UX, ...,
                sigma=sigma,
                at=at, leaveoneout=leaveoneout,
                weights=weights,
                edge=edge, diggle=diggle, positive=TRUE)
  #' extract smoothing bandwidth actually used
  sigma <- attr(LX, "sigma")
  varcov <- attr(LX, "varcov")
  #' initialise result
  Z <- LX
  Z[] <- NA
  #' guard against underflow
  tinythresh <- 8 * .Machine$double.eps
  if(underflow <- (min(LX) < tinythresh)) {
    Bad <- (LX < tinythresh)
    warning(paste("Numerical underflow detected at",
                  percentage(Bad, 1), "of", paste0(atName, ";"),
                  "sigma is probably too small"),
            call.=FALSE)
    #' apply l'Hopital's Rule at the problem locations
    Z[Bad] <- nnmark(X, at=at, xy=LX)[Bad]
    Good <- !Bad
  }
  #' compute
  for(k in 1:Num) {
    #' cumulative spatial weight of points with marks <= m_[k]
    if(k == Num) {
      Acum.k <- 1
    } else {
      w.k <- (m <= um[k]) * (weights %orifnull% 1)
      Lcum.k <- density(UX,
                        weights=w.k,
                        sigma=sigma, varcov=varcov, 
                        xy=LX,
                        at=at, leaveoneout=leaveoneout,
                        edge=edge, diggle=diggle, positive=TRUE)
      Acum.k <- Lcum.k/LX
    }
    if(k == 1) {
      #' region where quantile is um[1]
      relevant <- (Acum.k >= prob)
      if(underflow) relevant <- relevant & Good
      if(any(relevant)) {
        Z[relevant] <- um[1]
        if(verbose)
          splat("value um[1] =", um[1], "assigned to", sum(relevant), atName)
      }
    } else {
      #' region where quantile is between um[k-1] and um[k]
      unassigned <- (Acum.kprev < prob)
      if(underflow) unassigned <- unassigned & Good
      if(!any(unassigned)) break
      relevant <- unassigned & (Acum.k >= prob)
      if(any(relevant)) {
        if(type == 1) {
          ## left-continuous inverse
          left <- (Acum.k > prob)
          Z[relevant & left] <- um[k-1]
          Z[relevant & !left] <- um[k]
        } else if(type == 4) {
          ## linear interpolation
          Z[relevant] <- um[k-1] +
            (um[k] - um[k-1]) * ((prob - Acum.kprev)/(Acum.k - Acum.kprev))[relevant]
        }
        if(verbose)
          splat("values between",
                paste0("um", paren(k-1, "[")), "=", um[k-1],
                "and",
                paste0("um", paren(k, "[")), "=", um[k],
                "assigned to", sum(relevant), atName)
      }
    }
    Acum.kprev <- Acum.k
  }
  attr(Z, "sigma") <- sigma
  attr(Z, "varcov") <- varcov
  return(Z)
}
