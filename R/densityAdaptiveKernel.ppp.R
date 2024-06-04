#'
#'   densityAdaptiveKernel.ppp.R
#'
#'   $Revision: 1.16 $  $Date: 2024/06/04 03:09:11 $
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
  if(missing(bw) || is.null(bw)) {
    bw <- do.call.matched(bw.abram,
                          resolve.defaults(list(X=quote(X), at="points"),
                                           list(...)),
                          extrargs=names(args(as.mask)))
  } else if(is.numeric(bw)) {
    check.nvector(bw, nX, oneok=TRUE, vname="bw")
    if(length(bw) == 1) bw <- rep(bw, nX)
  } else if(is.im(bw)) {
    bw <- safelookup(bw, X, warn=FALSE)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of image 'bw'",
           call.=FALSE)
  } else if(inherits(bw, "funxy")) {
    bw <- bw(X)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of function 'bw'",
           call.=FALSE)
  } else stop("Argument 'bw' should be a numeric vector or a pixel image")

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
  return(ZZ)
}


densityAdaptiveKernel.ppplist <- 
densityAdaptiveKernel.splitppp <- function(X, bw=NULL, ...,
                                           weights=NULL) {
  n <- length(X)
  bw      <- ensure.listarg(bw,
                            n=n,
                            singletypes=c("NULL", "im", "funxy"),
                            xtitle="bw")
  weights <- ensure.listarg(weights,
                            n=n,
                            singletypes=c("NULL", "im", "funxy", "expression"),
                            xtitle="weights")
  y <- mapply(densityAdaptiveKernel.ppp, X=X, bw=bw, weights=weights,
              MoreArgs=list(...),
              SIMPLIFY=FALSE)
  return(as.solist(y, demote=TRUE))
}


## move this to spatstat.data when stable

ensure.listarg <- function(x, n, singletypes=character(0), 
                           xtitle=NULL, things="point patterns") {
  if(inherits(x, singletypes)) {
    ## single object: replicate it
    x <- rep(list(x), n)
    return(x)
  } 
  if(!is.list(x)) {
    ## error 
    if(is.null(xtitle)) xtitle <- short.deparse(substitute(x))
    whinge <- paste(xtitle, "should be a list")
    if(length(singletypes)) {
      otypes <- setdiff(singletypes, "NULL")
      if(length(otypes))
        whinge <- paste(whinge,
                        "or an object of class",
                        commasep(dQuote(otypes), "or"))
      if("NULL" %in% singletypes)
        whinge <- paste(whinge, "or NULL")
    }
    stop(whinge, call.=FALSE)
  }
  nx <- length(x)
  if(nx != n) {
    if(is.null(xtitle)) xtitle <- short.deparse(substitute(x))
    whinge <- paste("The length of",
                    sQuote(xtitle), 
                    "should equal the number of",
                    things,
                    paren(paste(nx, "!=", n)))
    stop(whinge, call.=FALSE)
  }
  return(x)
}
  
