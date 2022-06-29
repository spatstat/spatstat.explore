#'
#'   adaptive.density.R
#'
#'   $Revision: 1.3 $  $Date: 2022/06/29 03:05:22 $
#'

adaptive.density <- function(X, ..., method=c("voronoi", "kernel", "nearest")) {
  method <- match.arg(method)
  result <- switch(method,
                   voronoi = densityVoronoi(X, ...),
                   kernel  = densityAdaptiveKernel(X, ...),
                   nearest = {
                     if(is.lpp(X)) stop("not implemented for lpp objects")
                     nndensity(X, ...)
                   })
  return(result)
}
