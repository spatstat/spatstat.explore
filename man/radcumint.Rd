\name{radcumint}
\alias{radcumint}
\title{
  Radial Cumulative Integral
}
\description{
  Compute the cumulative integral of an image over increasing radial
  distances from the origin.
}
\usage{
radcumint(X, \dots, origin, Xname, result = c("fv", "im"))
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"})
    with numerical or logical values.
  }
  \item{\dots}{
    Ignored.
  }
  \item{origin}{
    Optional. Origin about which the rotations should be performed.
    Either a numeric vector or a character string as described
    in the help for \code{\link[spatstat.geom]{shift.owin}}.
  }
  \item{Xname}{
    Optional name for \code{X} to be used in the
    function labels.
  }
  \item{result}{
    Character string specifying the kind of result required:
    either a function object or a pixel image.
  }
}
\details{
  This command computes, for each possible distance \eqn{r},
  the integral of the pixel values lying inside the disc of radius
  \eqn{r} centred at the origin. 

  If \code{result="fv"} (the default) the result is a function
  object \code{f} of class \code{"fv"}. For each value of radius \eqn{r},
  the function value \code{f(r)} is the integral of \code{X}
  over the disc of radius \eqn{r}.

  If \code{result="im"} the result is a pixel image, with the same
  dimensions as \code{X}. At a given pixel, the result is
  equal to \code{f(r)} where \code{r} is the distance from the given
  pixel to the origin. That is, at any given pixel, the resulting value
  is the integral of \code{X} over the disc
  centred at the origin whose boundary passes through the given pixel. 
}
\value{
  An object of class \code{"fv"} or \code{"im"},
  with the same coordinate units as \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{rotmean}}, \code{\link{spatialcdf}}
}
\examples{
   D <- density(redwood)
   plot(radcumint(D))
   plot(radcumint(D, result="im"))
}
\keyword{spatial}
\keyword{math}
