\name{bw.bdh}
\alias{bw.bdh}
\title{
  Adjust Bandwidth Selection Rule to Account for Inhomogeneity
}
\description{
  Calculate a bandwidth for kernel estimation of the inhomogeneous pair
  correlation function. Given a bandwidth selection rule or bandwidth
  value which would be appropriate for a homogeneous point pattern,
  this function adjusts the bandwidth to account for inhomogeneity.
}
\usage{
bw.bdh(X, lambda=NULL, \dots, base=bw.stoyan, k=2)
}

\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function of \code{X}.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
    If \code{lambda} is missing or \code{NULL}, the intensity will
    be estimated from \code{X} using
    \code{\link[spatstat.explore]{density.ppp}}.
  }
  \item{\dots}{
    Arguments passed to
    \code{\link[spatstat.explore]{density.ppp}} to control the
    estimation of intensity (if \code{lambda=NULL}),
    or arguments passed to \code{base} if it is a function.
  }
  \item{base}{
    Bandwidth selection rule, or bandwidth value, that will be
    adjusted. Either a single numeric value, or a function that
    will be applied to the pattern \code{X} to obtain a bandwidth value,
    or a character string giving the name of such a function.
  }
  \item{k}{
    Integer exponent for calculating the adjustment factor.
  }
}
\details{
  This function calculates a bandwidth value for kernel estimation
  of the inhomogeneous pair correlation function
  using \code{\link{pcfinhom}}.

  It takes a bandwidth value or bandwidth selection rule
  that would be appropriate if the point process were homogeneous,
  adjusts the bandwidth to account for inhomogeneity,
  and returns the adjusted bandwidth value.

  The adjusted (inhomogeneous) bandwidth is the original (homogeneous) bandwidth
  multiplied by the \smoothpcfpapercite variance controlling adjustment factor.
  
  First a numerical bandwidth value, appropriate for a homogeneous process,
  is calculated. The default is to apply the extrapolated
  Stoyan rule-of-thumb \code{\link{bw.stoyan}} to the point pattern
  \code{X}. If \code{base} is specified, it may be either a 
  numeric value for the bandwidth, or another function that will be
  applied to \code{X} to calculate a bandwidth value.

  Next the intensity value at each point of \code{X} is evaluated.
  The argument \code{lambda} may be:
  \itemize{
    \item a numeric vector
    giving the intensity values at the points of \code{X};
    \item 
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations;
    \item
    a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
    The intensity of the fitted model will be evaluated at each
    point of \code{X}. By default, the fitted model is updated by
    re-fitting it to \code{X} before the intensity is evaluated.
    Updating can be disabled by setting \code{update=FALSE}.
    \item
    a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location;
    \item
    missing or \code{NULL}. In this case, the intensity will
    be estimated from \code{X} using 
    \code{\link[spatstat.explore]{density.ppp}}.
    Arguments \code{\dots} controlling the kernel estimation include 
    \code{sigma}, \code{varcov} and \code{kernel}.
  }
  
  Finally the bandwidth is adjusted by multiplying it by the
  \smoothpcfpapercite variance-controlling factor 
  \deqn{
    a = (n^{-1} \sum_i \lambda_i ) \, (n^{-1} \sum_i \lambda_i^{-k} )^{1/k}
  }{
    a = ((1/n) * sum[i] lambda[i]) ((1/n) * sum[i] (1/lambda[i]^k))^(1/k)
  }
  where \eqn{\lambda_i}{lambda[i]} is the value of \code{lambda}
  for the \eqn{i}th data point \code{X[i]}.

  When \code{k=2} (the default), the adjustment factor is
  \deqn{
    a = (n^{-1} \sum_i \lambda_i ) \, \sqrt{n^{-1} \sum_i \lambda_i^{-2}} )
  }{
    a = ((1/n) * sum[i] lambda[i]) sqrt((1/n) * sum[i] (1/lambda[i]^2))
  }
}
\value{
  A finite positive numerical value giving the selected bandwidth (the standard
  deviation of the smoothing kernel).
  The result has an attribute \code{"adjust"} giving the adjustment
  factor \eqn{a}.
}
\seealso{
  \code{\link{bw.stoyan}},
  \code{\link{bw.pcf}},
  \code{\link{pcfinhom}}
}
\examples{
  if(require(spatstat.model)) {
    fit <- ppm(japanesepines ~ x)
    (b <- bw.bdh(japanesepines, fit))
    attr(b, "adjust")
  }
}
\references{
  \smoothpcfpaper
}
\author{\adrian, \tilman and \martinH.}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
