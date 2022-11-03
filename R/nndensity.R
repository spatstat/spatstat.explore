#
#  nndensity.R
#
#  Density estimation based on nn distance
#
#  $Revision: 1.4 $  $Date: 2022/11/03 11:08:33 $
#

nndensity <- function(x, ...) {
  UseMethod("nndensity")
}

nndensity.ppp <- function(x, k, ..., verbose=TRUE) {
  if(missing(k) || is.null(k)) {
    k <- round(sqrt(npoints(x)))
    if(verbose) cat(paste("k=", k, "\n"))
  } else {
    check.1.integer(k)
    if(verbose && k == 1) warning("k=1 will produce strange results")
  }
  # distance to k-th nearest neighbour
  D <- nnmap(x, k=k, what="dist", ...)
  # area searched
  A <- eval.im(pi * D^2)
  # distance to boundary
  B <- bdist.pixels(as.owin(D))
  # handle edge effects
  edge <- solutionset(B < D)
  # centres of all pixels where edge effect occurs
  xy <- rasterxy.mask(edge, drop=TRUE)
  # corresponding values of distance
  rr <- D[edge, drop=TRUE]
  # compute actual search area
  X <- as.ppp(xy, W=as.owin(x), check=FALSE)
  A[edge] <- discpartarea(X, matrix(rr, ncol=1))
  # finally compute intensity estimate
  L <- eval.im(k/A)
  return(L)
}
