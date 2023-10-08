\name{adaptive.density}
\alias{adaptive.density}
\title{Adaptive Estimate of Intensity of Point Pattern}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern.
}
\usage{
adaptive.density(X, \dots, method=c("voronoi","kernel", "nearest"))
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"} or
    \code{"lpp"}).}
  \item{method}{Character string specifying the estimation method}
  \item{\dots}{
    Additional arguments passed to
    \code{\link{densityVoronoi}}, \code{\link{densityAdaptiveKernel}}
    or \code{\link{nndensity}}.
  }
}
\details{
  This function is an alternative to \code{\link{density.ppp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset. The result is a pixel image giving the estimated intensity.

  If \code{method="voronoi"} the data are passed to the function
  \code{\link{densityVoronoi}} which estimates the intensity using
  the Voronoi-Dirichlet tessellation.

  If \code{method="kernel"} the data are passed to the function
  \code{\link{densityAdaptiveKernel}} which estimates the intensity
  using a variable-bandwidth kernel estimator.

  If \code{method="nearest"} the data are passed to the function
  \code{\link{nndensity}} which estimates the intensity using the
  distance to the \code{k}-th nearest data point. (This is not supported
  when \code{X} has class \code{"lpp"}.)
}
\value{
  A pixel image (object of class \code{"im"}) whose values are
  estimates of the intensity of \code{X}.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{densityVoronoi}},
  \code{\link{densityAdaptiveKernel}},
  \code{\link{nndensity}},
  \code{\link{im.object}}.
}
\examples{
  plot(adaptive.density(nztrees, 1), main="Voronoi estimate")
}
\author{
  \spatstatAuthors and \mehdi.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Adaptive smoothing} 
