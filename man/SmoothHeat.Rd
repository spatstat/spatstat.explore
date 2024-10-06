\name{SmoothHeat}
\alias{SmoothHeat}
\title{Spatial Smoothing of Data by Diffusion}
\description{
  Generic function to perform spatial smoothing of spatial data
  by diffusion.
}
\usage{
  SmoothHeat(X, \dots)
}
\arguments{
  \item{X}{Some kind of spatial data}
  \item{\dots}{Arguments passed to methods.}
}
\details{
  This generic function calls an appropriate method
  to perform spatial smoothing on the spatial dataset \code{X}
  using diffusion.

  Methods for this function include
  \itemize{
    \item \code{\link[spatstat.explore]{SmoothHeat.ppp}} for point patterns
    \item \code{\link[spatstat.explore]{SmoothHeat.im}} for pixel images.
  }
}
\seealso{
  \code{\link[spatstat.explore]{SmoothHeat.ppp}},
  \code{\link[spatstat.explore]{SmoothHeat.im}}.
}
\value{
  An object containing smoothed values of the input data,
  in an appropriate format. See the documentation for the methods.
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
