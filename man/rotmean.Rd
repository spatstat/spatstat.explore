\name{rotmean}
\alias{rotmean}
\title{
  Rotational Average of a Pixel Image
}
\description{
  Compute the average pixel value over all rotations of the image
  about the origin, as a function of distance from the origin.
}
\usage{
rotmean(X, ..., origin, padzero=TRUE, Xname, result=c("fv", "im"), adjust=1)
}
\arguments{
  \item{X}{
    A pixel image.
  }
  \item{\dots}{
    Ignored.
  }
  \item{origin}{
    Optional. Origin about which the rotations should be performed.
    Either a numeric vector or a character string as described
    in the help for \code{\link[spatstat.geom]{shift.owin}}.
  }
  \item{padzero}{
    Logical. If \code{TRUE} (the default), the value of \code{X}
    is assumed to be zero outside the window of \code{X}.
    If \code{FALSE}, the value of \code{X} is taken to be undefined
    outside the window of \code{X}.
  }
  \item{Xname}{
    Optional name for \code{X} to be used in the
    function labels.
  }
  \item{result}{
    Character string specifying the kind of result required:
    either a function object or a pixel image.
  }
  \item{adjust}{
    Adjustment factor for bandwidth used in kernel smoothing.
  }
}
\details{
  This command computes, for each possible distance \eqn{r},
  the average pixel value of the pixels lying at
  distance \eqn{r} from the origin. Kernel smoothing is used
  to obtain a smooth function of \eqn{r}.

  If \code{result="fv"} (the default) the result is a function
  object of class \code{"fv"} giving the mean pixel value of \code{X}
  as a function of distance from the origin.

  If \code{result="im"} the result is a pixel image, with the same
  dimensions as \code{X}, giving the mean value of \code{X}
  over all pixels lying at the same distance from the origin
  as the current pixel.

  If \code{padzero=TRUE} (the default), the value of \code{X}
  is assumed to be zero outside the window of \code{X}. The rotational
  mean at a given distance \eqn{r} is the average value of the image
  \code{X} over the \emph{entire} circle of radius \eqn{r},
  including zero values outside the window if the circle
  lies partly outside the window.
  
  If \code{padzero=FALSE}, the value of \code{X} is taken to be
  undefined outside the window of \code{X}. The rotational mean
  is the average of the \code{X} values over the \emph{subset} of the circle
  of radius \eqn{r} that lies entirely inside the window.
}
\value{
  An object of class \code{"fv"} or \code{"im"},
  with the same coordinate units as \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{radcumint}}
}
\examples{
  online <- interactive()
  resolution <- if(online) 128 else 32
  Z <- setcov(square(1), dimyx=resolution)
  f <- rotmean(Z)
  if(online) {
    plot(rotmean(Z))
    plot(rotmean(Z, result="im"))
  }
}
\keyword{spatial}
\keyword{math}
