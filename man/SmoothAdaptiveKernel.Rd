\name{SmoothAdaptiveKernel}
\alias{SmoothAdaptiveKernel}
\alias{SmoothAdaptiveKernel.ppp}
\title{Spatial Smoothing of Marks using Adaptive Kernel}
\description{
  Performs spatial smoothing of numeric values observed
  at a set of irregular locations, using adaptive (variable-bandwidth)
  kernel smoothing.
}
\usage{
SmoothAdaptiveKernel(X, bw, \dots)

\method{SmoothAdaptiveKernel}{ppp}(X, bw, \dots,
       weights=NULL,       
       at=c("pixels", "points"))
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{bw}{
    Numeric vector of smoothing bandwidths for each point in \code{X},
    or a pixel image giving the smoothing bandwidth at each spatial
    location, or a spatial function of class \code{"funxy"} giving the
    smoothing bandwidth at each location.
    The default is to compute bandwidths using \code{\link{bw.abram.ppp}}.
  }
  \item{\dots}{
    Arguments passed to
    \code{\link{bw.abram.ppp}} to compute the smoothing bandwidths
    if \code{bw} is missing,
    or passed to \code{\link[spatstat.geom]{as.mask}}
    to control the spatial resolution of the result,
    or passed to \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}} to control
    calculations such as edge corrections.
  }
  \item{weights}{
    Optional vector of numeric weights for the points of \code{X}.
  }
  \item{at}{
    String specifying whether to compute the smoothed values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
}
\details{
  This function performs spatial smoothing of the
  mark values of \code{X} using spatially-adaptive kernel smoothing.

  The function \code{SmoothAdaptiveKernel} is generic.
  This file documents the method for point patterns,
  \code{SmoothAdaptiveKernel.ppp}.

  The argument \code{X} should be a marked point pattern (an object
  of class \code{"ppp"} with mark values given by \code{marks(X)}).
  The mark values will be interpreted as numerical observations
  recorded at the spatial locations of the points of \code{X}.
  This function applies spatial smoothing to the mark values.
  
  The argument \code{bw} specifies the smoothing bandwidths to be
  applied to each of the points in \code{X}. It may be a numeric vector
  of bandwidth values, or a pixel image or function yielding the
  bandwidth values. See \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}} for more
  information.

  The smoother is an adaptive (variable-bandwidth) version of the
  Nadaraya-Watson smoother.
  If the points of \code{X} are \eqn{x_1,\ldots,x_n}{x[1], ..., x[n]},
  the mark values are \eqn{v_1,\ldots,v_n}{v[1], ..., v[n]}
  and the corresponding bandwidths are
  \eqn{\sigma_1,\ldots,\sigma_n}{\sigma[1], ..., \sigma[n]}
  then the adaptive kernel smoother at a location \eqn{u} is
  \deqn{
    S(u) = \frac{
      \sum_{i=1}^n v_i k(u, x_i, \sigma_i)
    }{
      \sum_{i=1}^n k(u, x_i, \sigma_i)
    }
  }{
    S(u) = (sum[i] v[i] * e(x[i], k(u, x[i], \sigma[i])))/
           (sum[i] e(x[i], k(u, x[i], \sigma[i])))
  }
  where \eqn{k(u, v, \sigma)} is the value at \eqn{u}
  of the (possibly edge-corrected) smoothing kernel with bandwidth \eqn{\sigma}
  induced by a data point at \eqn{v}.

  The algorithm uses \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}}
  to compute the numerator and denominator of this expression.
}
\value{
  If \code{at="pixels"} (the default), the result is a pixel image.
  If \code{at="points"}, the result is a numeric vector with one entry
  for each data point in \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{bw.abram.ppp}},
  \code{\link[spatstat.explore]{Smooth.ppp}},
  \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}}.
}
\examples{
  X <- redwood3
  marks(X) <- X$x
  Z <- SmoothAdaptiveKernel(X, h0=0.1)
  plot(Z, main="Adaptive smoother")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Adaptive smoothing}
