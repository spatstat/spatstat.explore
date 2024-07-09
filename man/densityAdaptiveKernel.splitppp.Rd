\name{densityAdaptiveKernel.splitppp}
\alias{densityAdaptiveKernel.splitppp}
\alias{densityAdaptiveKernel.ppplist}
\title{Adaptive Kernel Estimate of Intensity for Split Point Pattern}
\description{
  Computes an adaptive estimate of the intensity function 
  (using a variable-bandwidth smoothing kernel)
  for each of the components of a split point pattern,
  or each of the point patterns in a list.
}
\usage{
\method{densityAdaptiveKernel}{splitppp}(X, bw=NULL, \dots, weights=NULL)

\method{densityAdaptiveKernel}{ppplist}(X, bw=NULL, \dots, weights=NULL)
}
\arguments{
  \item{X}{
    Split point pattern (object of class \code{"splitppp"}
    created by \code{\link[spatstat.geom]{split.ppp}}) to be smoothed.
    Alternatively a list of point patterns,
    of class \code{"ppplist"}.
  }
  \item{bw}{
    Smoothing bandwidths. See Details.
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link{densityAdaptiveKernel.ppp}}.
    These may include arguments that will be passed to
    \code{\link{bw.abram.ppp}} to compute
    the smoothing bandwidths if \code{bw} is missing,
    and arguments passed to \code{\link[spatstat.geom]{as.mask}}
    to control the spatial resolution of the result.
  }
  \item{weights}{
    Numerical weights for the points. See Details.
  }
}
\details{
  This function computes a spatially-adaptive kernel estimate of the
  spatially-varying intensity for each of the point patterns
  in the list \code{X}, using \code{\link{densityAdaptiveKernel.ppp}}.

  The argument \code{bw} specifies smoothing bandwidths
  for the data points.
  Normally it should be a list, with the same length as
  \code{x}. The entry \code{bw[[i]]} will determine the
  smoothing bandwidths for the pattern \code{x[[i]]}, and may be given in
  any format acceptable to \code{\link{densityAdaptiveKernel.ppp}}.
  For example, \code{bw[[i]]} can be
  a numeric vector of length equal to \code{npoints(x[[i]])},
  a single numeric value, 
  a pixel image (object of class \code{"im"}),
  an \code{expression}, or a function of class \code{"funxy"}.
  For convenience, \code{bw} can also be a single \code{expression},
  or a single pixel image, or a single function.
  If \code{bw} is missing or \code{NULL}, the default is to compute
  bandwidths using \code{\link{bw.abram.ppp}}.

  The argument \code{weights} specifies numerical case weights
  for the data points.
  Normally it should be a list, with the same length as
  \code{x}. The entry \code{weights[[i]]} will determine the
  case weights for the pattern \code{x[[i]]}, and may be given in
  any format acceptable to \code{\link{density.ppp}}.
  For example, \code{weights[[i]]} can be
  a numeric vector of length equal to \code{npoints(x[[i]])},
  a single numeric value, a numeric matrix,
  a pixel image (object of class \code{"im"}),
  an \code{expression}, or a function of class \code{"funxy"}.
  For convenience, \code{weights} can also be a single \code{expression},
  or a single pixel image (object of class \code{"im"}),
  or a single function of class \code{"funxy"}.
  If \code{weights} is missing or \code{NULL}, all weights are assumed
  to be equal to 1.
}
\value{
  A list of pixel images (objects of class \code{"im"})
  which can be plotted or printed;
  or a list of numeric vectors giving the values at specified points.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{densityAdaptiveKernel.ppp}},
  \code{\link{bw.abram.ppp}}.
}
\examples{
  X <- amacrine
  if(!interactive()) X <- X[c(TRUE,FALSE,FALSE,FALSE)]
  Z <- densityAdaptiveKernel(split(X), h0=0.15)
  plot(Z, main="Adaptive kernel estimate")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Adaptive smoothing}
