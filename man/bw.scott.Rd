\name{bw.scott}
\alias{bw.scott}
\alias{bw.scott.iso}
\title{
  Scott's Rule for Bandwidth Selection for Kernel Density
}
\description{
  Use Scott's rule of thumb to determine the smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.scott(X, isotropic=FALSE, d=NULL)

   bw.scott.iso(X)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"},
    \code{"lpp"}, \code{"pp3"} or \code{"ppx"}).
  }
  \item{isotropic}{
    Logical value indicating whether to compute a single
    bandwidth for an isotropic Gaussian kernel (\code{isotropic=TRUE})
    or separate bandwidths for each coordinate axis
    (\code{isotropic=FALSE}, the default).
  }
  \item{d}{
    Advanced use only.
    An integer value that should be used in Scott's formula
    instead of the true number of spatial dimensions.
  }
}
\details{
  These functions select a bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}
  or other appropriate functions.
  They can be applied to a point pattern
  belonging to any class \code{"ppp"}, \code{"lpp"}, \code{"pp3"}
  or \code{"ppx"}.

  The bandwidth \eqn{\sigma}{\sigma} is computed by the rule of thumb
  of Scott (1992, page 152, equation 6.42). 
  The bandwidth is proportional to \eqn{n^{-1/(d+4)}}{n^(-1/(d+4))}
  where \eqn{n} is the number of points and \eqn{d} is the number of
  spatial dimensions.

  This rule is very fast to compute. It typically produces a larger bandwidth
  than \code{\link{bw.diggle}}. It is useful for estimating
  gradual trend.

  If \code{isotropic=FALSE} (the default), \code{bw.scott} provides a
  separate bandwidth for each coordinate axis, and the result of the
  function is a vector, of length equal to the number of coordinates.
  If \code{isotropic=TRUE}, a single bandwidth value is computed
  and the result is a single numeric value.

  \code{bw.scott.iso(X)} is equivalent to
  \code{bw.scott(X, isotropic=TRUE)}.

  The default value of \eqn{d} is as follows:
  \tabular{ll}{
    \bold{class} \tab \bold{dimension} \cr
    \code{"ppp"} \tab 2 \cr
    \code{"lpp"} \tab 1 \cr
    \code{"pp3"} \tab 3 \cr
    \code{"ppx"} \tab number of spatial coordinates
  }
  The use of \code{d=1} for point patterns on a linear network
  (class \code{"lpp"}) was proposed by McSwiggan et al (2016)
  and Rakshit et al (2019).
}
\value{
  A numerical value giving the selected
  bandwidth, or a numerical vector giving the
  selected bandwidths for each coordinate.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.diggle}},
  \code{\link{bw.ppl}},
  \code{\link{bw.CvL}},
  \code{\link{bw.frac}}.
}
\examples{
  hickory <- split(lansing)[["hickory"]]
  b <- bw.scott(hickory)
  b
  if(interactive()) {
   plot(density(hickory, b))
  }
  bw.scott.iso(hickory)
  bw.scott(osteo$pts[[1]])
}
\references{
  Scott, D.W. (1992)
  \emph{Multivariate Density Estimation. Theory, Practice and
    Visualization}. 
  New York: Wiley.

}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Bandwidth selection}
