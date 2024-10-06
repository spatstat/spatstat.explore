\name{SmoothHeat.ppp}
\alias{SmoothHeat.ppp}
\title{Spatial Smoothing of Observations using Diffusion Estimate of Density}
\description{
  Performs spatial smoothing of numeric values observed
  at a set of irregular locations, using the diffusion estimate
  of the density.
}
\usage{
\method{SmoothHeat}{ppp}(X, sigma, \dots, weights=NULL)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
    with numeric marks.
  }
  \item{sigma}{
    Smoothing bandwidth. A single number giving the equivalent
    standard deviation of the smoother.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{densityHeat}}
    controlling the estimation of each marginal intensity,
    or passed to \code{\link[spatstat.geom]{pixellate.ppp}}
    controlling the pixel resolution.
  }
  \item{weights}{Optional numeric vector of weights associated with each
    data point.
  }
}
\details{
  This is the analogue of the Nadaraya-Watson smoother, using the
  diffusion smoothing estimation procedure (Baddeley et al, 2022).
  The numerator and denominator of the Nadaraya-Watson smoother are
  calculated using \code{\link[spatstat.explore]{densityHeat.ppp}}.
}
\value{
  Pixel image (object of class \code{"im"}) giving the smoothed
  mark value.
}
\seealso{
  \code{\link[spatstat.explore]{Smooth.ppp}} for the usual kernel-based
  smoother (the Nadaraya-Watson smoother)
  and \code{\link[spatstat.explore]{densityHeat}} for the diffusion estimate of density.
}
\author{
  \adrian, \tilman and Suman Rakshit.
}
\examples{
   plot(SmoothHeat(longleaf, 10))
}
\references{
  Baddeley, A., Davies, T., Rakshit, S., Nair, G. and McSwiggan, G. (2022)
  Diffusion smoothing for spatial point patterns.
  \emph{Statistical Science} \bold{37}, 123--142.
}
\keyword{spatial}
\keyword{smooth}
