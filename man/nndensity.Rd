\name{nndensity.ppp}
\alias{nndensity}
\alias{nndensity.ppp}
\title{
  Estimate Intensity of Point Pattern Using Nearest Neighbour Distances
}
\description{
  Estimates the intensity of a point pattern
  using the distance from each spatial location
  to the \code{k}th nearest data point.
}
\usage{
nndensity(x, ...)

\method{nndensity}{ppp}(x, k, ..., verbose = TRUE)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"})
    or some other spatial object.
  }
  \item{k}{
    Integer. The distance to the \code{k}th nearest data point
    will be computed. There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{nnmap}} and \code{\link[spatstat.geom]{as.mask}}
    controlling the pixel resolution.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, print the value of \code{k}
    when it is automatically selected. If \code{FALSE}, remain silent.
  }
}
\details{
  This function computes a quick estimate of the intensity of the point
  process that generated the point pattern \code{x}.

  For each spatial location \eqn{s}, let \eqn{d(s)} be the distance from \eqn{s}
  to the \eqn{k}-th nearest point in the dataset \code{x}.
  If the data came from a homogeneous
  Poisson process with intensity \eqn{\lambda}{lambda},
  then \eqn{\pi d(s)^2}{pi * d(s)^2} would follow a 
  negative exponential distribution with mean
  \eqn{1/\lambda}{1/lambda}, and the maximum likelihood estimate of
  \eqn{\lambda}{lambda} would be \eqn{1/(\pi d(s)^2)}{1/(pi * d(s)^2)}.
  This is the estimate computed by \code{nndensity},
  apart from an edge effect correction.

  See Cressie (1991, equation (8.5.14), p. 654)
  and Silverman (1986, p. 96).

  This estimator of intensity is relatively fast to compute, 
  and is spatially adaptive
  (so that it can handle wide variation in the intensity
  function). However, it implicitly assumes the points are independent,
  so it does not perform well if the pattern is strongly clustered
  or strongly inhibited.

  In normal use, the value of \code{k} should be at least 3.
  (Theoretically the estimator has infinite expected value if \eqn{k=1},
  and infinite variance if \eqn{k=2}.
  The computed intensity estimate will have infinite peaks
  around each data point if \code{k = 1}.)
  The default value of \code{k} is the square root of the number of
  points in \code{x}, which seems to work well in many cases.
  
  The window of \code{x} is digitised using \code{\link[spatstat.geom]{as.mask}}
  and the values \eqn{d(s)} are computed using \code{\link[spatstat.geom]{nnmap}}.
  To control the pixel resolution, see \code{\link[spatstat.geom]{as.mask}}.
}
\value{
  A pixel image (object of class \code{"im"}) giving the
  estimated intensity of the point process at each spatial location.
  Pixel values are intensities (number of points per unit area).
}
\references{
  Cressie, N.A.C. (1991)
  \emph{Statistics for spatial data}.
  John Wiley and Sons, New York.

  Silverman, B.W. (1986)
  \emph{Density Estimation}.
  Chapman and Hall, New York.
}
\seealso{
 \code{\link{density.ppp}},
 \code{\link[spatstat.geom]{intensity.ppp}}
 for alternative estimates of point process intensity.
}
\examples{
   plot(nndensity(swedishpines))
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
