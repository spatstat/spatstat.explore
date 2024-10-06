\name{blurHeat}
\alias{blurHeat}
\alias{blurHeat.im}
\alias{SmoothHeat.im}
\title{
  Diffusion Blur
}
\description{
  Blur a Pixel Image by Applying Diffusion
}
\usage{
blurHeat(X, \dots)

\method{blurHeat}{im}(X, sigma, \dots,
     connect = 8, symmetric = FALSE, k= 1, show = FALSE)

\method{SmoothHeat}{im}(X, sigma, \dots)

}
\arguments{
  \item{X}{
    Pixel image (object of class \code{"im"}).
  }
  \item{sigma}{
    Smoothing bandwidth. A numeric value, a pixel image or
    a \code{function(x,y)}.
  }
  \item{\dots}{
    Ignored by \code{blurHeat.im}.
  }
  \item{connect}{
    Grid connectivity: either 4 or 8.
  }
  \item{symmetric}{
    Logical value indicating whether to \emph{force} the algorithm
    to use a symmetric random walk.
  }
  \item{k}{
    Integer. Calculations will be performed by repeatedly multiplying
    the current state by the \code{k}-step transition matrix.
  }
  \item{show}{
    Logical value indicating whether to plot successive iterations.
  }
}
\details{
  The function \code{blurHeat} is generic. 

  This help file documents the method \code{blurHeat.im} for pixel images
  (objects of class \code{"im"}). This is currently equivalent
  to \code{SmoothHeat.im}, which is also documented here.
  
  If \code{sigma} is a numeric value, then
  the classical time-dependent heat equation is solved
  up to time \code{t = sigma^2} starting with the initial
  condition given by the image \code{X}. This has the effect
  of blurring the input image \code{X}.

  If \code{sigma} is a function or a pixel image, then
  it is treated as a spatially-variable diffusion rate,
  and the corresponding heat equation is solved.

  This command can be used to calculate the expected value
  of the diffusion estimator of intensity (\code{\link[spatstat.explore]{densityHeat}})
  when the true intensity is known.
}
\value{
  A pixel image on the same raster as \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{densityHeat}},
  \code{\link[spatstat.explore]{blur}}.
}
\examples{
  Z <- as.im(function(x,y) { sin(10*x) + sin(9*y) }, letterR)
  ZZ <- blurHeat(Z, 0.2)
  plot(solist(original=Z, blurred=ZZ), main="")
}
\keyword{spatial}
\keyword{math}
