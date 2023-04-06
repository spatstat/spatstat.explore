##
## spatialcdf.R
##
##  $Revision: 1.10 $ $Date: 2023/04/06 00:14:21 $
##

spatialcdf <- function(Z, weights=NULL, normalise=FALSE, ...,
                       W=NULL, Zname=NULL) {
  Zdefaultname <- singlestring(short.deparse(substitute(Z)))
  if(is.character(Z) && length(Z) == 1) {
    if(is.null(Zname)) Zname <- Z
    switch(Zname,
           x={
             Z <- function(x,y) { x }
           }, 
           y={
             Z <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
  }
  if(is.null(Zname)) Zname <- Zdefaultname
  ##
  if(inherits(weights, c("ppm", "kppm", "dppm"))) {
    model <- weights
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required", call.=FALSE)
    df <- spatstat.model::spatialCovariateUnderModel(model, Z)
    G <- with(df, ewcdf(Z, wt, normalise=normalise))
    wtname <- if(normalise) "fraction of points" else "number of points"
  } else {
    if(is.null(W)) W <- as.owin(weights, fatal=FALSE)
    if(is.null(W)) W <- as.owin(Z, fatal=FALSE)
    if(is.null(W)) stop("No information specifying the spatial window")
    M <- as.mask(W, ...)
    loc <- as.ppp(rasterxy.mask(M, drop=TRUE), W=W, check=FALSE)
    pixelarea <- with(unclass(M), xstep * ystep)
    if(is.null(weights)) {
      Zvalues <- evaluateCovariateAtPoints(Z, loc, ...)
      G <- ewcdf(Zvalues, normalise=normalise, adjust=pixelarea)
      wtname <- if(normalise) "fraction of area" else "area"
    } else {
      Zvalues <- evaluateCovariateAtPoints(Z, loc, ...)
      wtvalues <- evaluateCovariateAtPoints(weights, loc, ...)
      G <- ewcdf(Zvalues, wtvalues, normalise=normalise, adjust=pixelarea)
      wtname <- if(normalise) "fraction of weight" else "weight"
    }
  }
  class(G) <- c("spatialcdf", class(G))
  attr(G, "call") <- sys.call()
  attr(G, "Zname") <- Zname
  attr(G, "ylab") <- paste("Cumulative", wtname)
  return(G)
}

plot.spatialcdf <- function(x, ..., xlab, ylab, do.points=FALSE) {
  if(missing(xlab) || is.null(xlab))
    xlab <- attr(x, "Zname")
  if(missing(ylab) || is.null(ylab))
    ylab <- attr(x, "ylab")
  if(inherits(x, "ecdf")) {
    plot.ecdf(x, ..., xlab=xlab, ylab=ylab, do.points=do.points)
  } else {
    plot.stepfun(x, ..., xlab=xlab, ylab=ylab, do.points=do.points)
  }
}

