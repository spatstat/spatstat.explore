#'
#'   densityAdaptiveKernel.ppp.R
#'
#'   $Revision: 1.19 $  $Date: 2026/09/19 07:50:50 $
#'
#'
#'  Adaptive kernel smoothing via 3D FFT
#'

densityAdaptiveKernel.ppp <- function(X, bw, ...,
                                      weights=NULL,
                                      at=c("pixels", "points"),
                                      edge=TRUE, 
                                      ngroups) {
  stopifnot(is.ppp(X))
  at <- match.arg(at)
  nX <- npoints(X)

  if("se" %in% names(list(...)))
   stop("Standard errors are not yet supported in densityAdaptiveKernel.ppp",
        call.=FALSE)
  
  if(nX == 0)
    switch(at,
           points = return(numeric(nX)),
           pixels = return(as.im(0, W=Window(X), ...)))
                     
  if(missing(ngroups) || is.null(ngroups)) {
    ## default rule
    ngroups <- max(1L, floor(sqrt(nX)))
  } else if(any(is.infinite(ngroups))) {
    ngroups <- nX
  } else {
    check.1.integer(ngroups)
    ngroups <- min(nX, ngroups)
  }

  if(weighted <- !is.null(weights)) {
    check.nvector(weights, nX, oneok=TRUE, vname="weights")
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else weights <- rep(1,nX)

  ## determine bandwidth for each data point
  if(missing(bw)) bw <- NULL
  bw <- resolve.adaptive.bandwidths(X, bw, ...)

  #' divide bandwidths into groups
  if(ngroups == nX) {
    ## every data point is a separate group
    groupid <- 1:nX
    qmid <- bw
  } else {
    ## usual case
    p <- seq(0,1,length=ngroups+1)
    qbands <- quantile(bw, p)
    groupid <- findInterval(bw,qbands,all.inside=TRUE)
    #' map to middle of group
    pmid <- (p[-1] + p[-length(p)])/2
    qmid   <- quantile(bw, pmid)
  }

  marks(X) <- if(weighted) weights else NULL
  group <- factor(groupid, levels=1:ngroups)
  Y <- split(X, group)

  Z <- mapply(density.ppp,
              x=Y,
              sigma=as.list(qmid),
              weights=lapply(Y, marks),
              MoreArgs=list(edge=edge, at=at, ...),
              SIMPLIFY=FALSE)

  ZZ <- switch(at,
               pixels = im.apply(Z, "sum"),
               points = unsplit(Z, group))
  attr(ZZ, "bw") <- bw
  return(ZZ)
}


densityAdaptiveKernel.ppplist <- 
densityAdaptiveKernel.splitppp <- function(X, bw=NULL, ...,
                                           weights=NULL) {
  n <- length(X)
  bw      <- ensure.nlist(bw,
                          n=n,
                          singletypes=c("NULL", "im", "funxy"),
                          xtitle="bw")
  weights <- ensure.nlist(weights,
                          n=n,
                          singletypes=c("NULL", "im", "funxy", "expression"),
                          xtitle="weights")
  y <- mapply(densityAdaptiveKernel.ppp, X=X, bw=bw, weights=weights,
              MoreArgs=list(...),
              SIMPLIFY=FALSE)
  return(as.solist(y, demote=TRUE))
}

resolve.adaptive.bandwidths <- function(X, bw=NULL, ...,
                                        adjust=1, sigma=NULL, warn=FALSE) {
  ## determine bandwidth for each data point of X
  if(!is.null(sigma)) {
    if(is.null(bw)) {
      ## catch inadvertent use of wrong argument name
      bw <- sigma
    } else {
      ## conflict
      if(warn)
        warning("Both arguments bw and sigma were specified; ignoring sigma",
                call.=FALSE)
    }
  }
  if(is.null(bw)) {
    ## no bandwidth information given
    bw <- do.call.matched(bw.abram,
                          resolve.defaults(list(X=quote(X), at="points"),
                                           list(...)),
                          extrargs=names(args(as.mask)))
  } else if(is.numeric(bw)) {
    ## bandwidth value for each data point
    nX <- npoints(X)
    check.nvector(bw, nX, oneok=TRUE, vname="bw")
    if(length(bw) == 1) bw <- rep(bw, nX)
  } else if(is.im(bw)) {
    ## bandwidth for a grid of spatial locations -- look up at data points
    bw <- safelookup(bw, X, warn=FALSE)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of image 'bw'",
           call.=FALSE)
  } else if(inherits(bw, "funxy")) {
    ## bandwidth for any spatial location --- evaluate at data points
    bw <- bw(X)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of function 'bw'",
           call.=FALSE)
  } else stop(paste("Argument 'bw' should be a numeric vector,",
                    "a pixel image, or a function(x,y)"),
              call.=FALSE)
  if(!missing(adjust)) {
    check.1.real(adjust)
    stopifnot(adjust > 0)
    bw <- adjust * bw
  }
  return(bw)
}

