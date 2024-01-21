\name{SpatialMedian.ppp}
\alias{SpatialMedian.ppp}
\title{
  Spatially Weighted Median of Values at Points
}
\description{
  Given a spatial point pattern with numeric marks,
  compute a weighted median of the mark values,
  with spatially-varying weights that depend on distance to the data points.
}
\usage{
 \method{SpatialMedian}{ppp}(X, sigma = NULL, \dots,
        type = 4, at = c("pixels", "points"), leaveoneout = TRUE,
        weights = NULL, edge = TRUE, diggle = FALSE, verbose = FALSE)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"})
    with numeric marks.
  }
  \item{sigma}{
    Smoothing bandwidth, passed to \code{\link{density.ppp}}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{density.ppp}} controlling the
    spatial smoothing.
  }
  \item{type}{
    Integer specifying the type of median
    (using the convention of \code{\link[stats]{quantile.default}};
    see Details).
    Only types 1 and 4 are currently implemented.
  }
  \item{at}{
    Character string indicating whether to compute the median
    at every pixel of a pixel image (\code{at="pixels"}, the default)
    or at every data point of \code{X} (\code{at="points"}).
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
  \item{weights}{
    Optional vector of numeric weights attached to the points of \code{X}.
  }
  \item{edge,diggle}{
    Arguments passed to \code{\link{density.ppp}} to
    determine the edge correction.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports
    during the calculation.
  }
}
\details{
  The argument \code{X} should be a spatial point pattern
  (object of class \code{"ppp"}) with numeric marks.

  The algorithm computes the weighted median of the mark values
  at each desired spatial location, using spatially-varying weights
  which depend on distance to the data points.

  Suppose the data points are at spatial locations
  \eqn{x_1,\ldots,x_n}{x[1], ..., x[n]}
  and have mark values 
  \eqn{y_1,\ldots,y_n}{y[1], ..., y[n]}.
  For a query location \eqn{u}, the smoothed median is defined
  as the weighted median of the mark values
  \eqn{y_1,\ldots,y_n}{y[1], ..., y[n]} with weights
  \eqn{w_1,\ldots,w_n}{w[1], ..., w[n]},
  where
  \deqn{
    w_i = \frac{k(u,x_i)}{\sum_{j=1}^n k(u,x_j)}
  }{
    w[i] = k(u,x[i])/(k(u, x[1]) + ... + k(u, x[n]))
  }
  where \eqn{k(u,v)} is the smoothing kernel with bandwidth \code{sigma}

  If \code{at="points"} and \code{leaveoneout=TRUE}, then
  a leave-one-out calculation is performed, which means that
  when the query location is a data point \eqn{x_i}{x[i]},
  the value at the data point is ignored, and 
  the weighted median is computed from the values \eqn{y_j}{y[j]}
  for all \eqn{j} not equal to \eqn{i}.
}
\value{
  \emph{If \code{X} has a single column of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is
    a pixel image (object of class \code{"im"}). 
    \item
    If \code{at="points"}, the result is a numeric vector
    of length equal to the number of points in \code{X}.
  }
  \emph{If \code{X} has a data frame of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is a named list of 
    pixel images (object of class \code{"im"}). There is one
    image for each column of marks. This list also belongs to
    the class \code{"solist"}, for which there is a plot method.
    \item
    If \code{at="points"}, the result is a data frame
    with one row for each point of \code{X},
    and one column for each column of marks. 
    Entries are values of the interpolated function at the points of \code{X}.
  }
  The return value has attributes
  \code{"sigma"} and \code{"varcov"} which report the smoothing
  bandwidth that was used.

  The calculation of the median value depends on the argument
  \code{type} which is interpreted in the same way as for
  \code{\link[stats]{quantile.default}}. Currently, only types 1 and 4
  are implemented. 
  If \code{type=1}, the median is always
  one of the mark values (one of the values in \code{marks(x)}).
  If \code{type=4} (the default),
  the median value is obtained by linearly interpolating
  between mark values.
  Note that the default values of \code{type}
  in \code{SpatialMedian.ppp} and \code{\link{SpatialQuantile.ppp}}
  are different.
}
\author{
  \adrian.
}
\seealso{
  Generic function \code{\link{SpatialMedian}}.

  \code{\link{SpatialQuantile}} and \code{\link{SpatialQuantile.ppp}}
  for other quantiles.
  
  \code{\link{Smooth.ppp}} for the spatially weighted average.
}
\examples{
   X <- longleaf
   if(!interactive()) {
    ## mark values rounded to nearest multiple of 10 to reduce check time
    marks(X) <- round(marks(X), -1)
   }
   Z <- SpatialMedian(X, sigma=30)
   ZX <- SpatialMedian(X, sigma=30, at="points")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
