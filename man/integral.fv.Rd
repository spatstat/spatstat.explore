\name{integral.fv}
\alias{integral.fv}
\title{
  Compute Integral of Function Object
}
\description{
  Compute the integral of a function 
  over a specified range.
}
\usage{
\method{integral}{fv}(f, domain = NULL, ...)
}
\arguments{
  \item{f}{
    A function value table
    (object of class \code{"fv"}).
  }
  \item{domain}{
    Optional. Range of values of the argument \eqn{x} over which the
    density \eqn{f(x)} should be integrated. A numeric vector of length
    2 giving the minimum and maximum values of \eqn{x}. Infinite limits
    are permitted.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link[spatstat.univar]{integral}}.
  It computes the numerical integral
  \deqn{
    I = \int f(x) dx
  }{
    I = integral f(x) dx
  }
  of the function object \code{f}.
  If \code{domain} is specified, the integral is restricted to the
  interval of \eqn{x} values given by the \code{domain}.

  The result is a numeric value or numeric vector containing one entry
  for each column of function values in \code{f}.
  
  Integrals are calculated numerically using the trapezoidal rule
  restricted to the domain given.
}
\value{
  A single numerical value, or a numeric vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{fv.object}}, \code{\link[spatstat.univar]{integral}}

  \code{\link[spatstat.univar]{stieltjes}}
}
\examples{
  g <- pcf(redwood, divisor="d")
  integral(g, domain=c(0, 0.1))
}
\keyword{methods}
\keyword{univar}
\keyword{nonparametric}
\keyword{math}
