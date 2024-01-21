\name{SpatialQuantile}
\alias{SpatialMedian}
\alias{SpatialQuantile}
\title{
  Spatially Weighted Median or Quantile
}
\description{
  Compute a weighted median or weighted quantile of spatial data.
}
\usage{
 SpatialMedian(X, \dots)

 SpatialQuantile(X, prob = 0.5, \dots)
}
\arguments{
  \item{X}{
    A spatial data object.
  }
  \item{prob}{
    Probability for which the quantile is required.
    A single numeric value between 0 and 1.
    Default is to calculate the median.
  }
  \item{\dots}{
    Further arguments passed to methods.
  }
}
\details{
  The functions \code{SpatialMedian} and \code{SpatialQuantile}
  are generic. They calculate spatially weighted medians
  and quantiles of spatial data. The details depend on the class of
  \code{X}.

  There are methods for spatial point patterns (class \code{"ppp"})
  and possibly for other objects.
}
\author{
  \adrian.
}
\seealso{
  Methods \code{\link{SpatialMedian.ppp}}, \code{\link{SpatialQuantile.ppp}}.
  
  \code{\link{Smooth}} for the spatially weighted average.
}
\keyword{spatial}
\keyword{smooth}
