\name{bw.smoothppp}
\alias{bw.smoothppp}
\title{
  Cross Validated Bandwidth Selection for Spatial Smoothing
}
\description{
  Uses least-squares cross-validation to select a smoothing bandwidth
  for spatial smoothing of marks.
}
\usage{
   bw.smoothppp(X, nh = spatstat.options("n.bandwidth"),
                hmin=NULL, hmax=NULL, warn=TRUE, kernel="gaussian",
                varcov1=NULL)
}

\arguments{
  \item{X}{
    A marked point pattern with numeric marks.
  }
  \item{nh}{
    Number of trial values of smoothing bandwith \code{sigma}
    to consider. The default is 32.
  }
  \item{hmin, hmax}{
    Optional. Numeric values.
    Range of trial values of smoothing bandwith \code{sigma}
    to consider. There is a sensible default.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the minimum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
  \item{kernel}{
    The smoothing kernel.
    A character string specifying the smoothing kernel
    (current options are \code{"gaussian"}, \code{"epanechnikov"},
    \code{"quartic"} or \code{"disc"}).
  }
  \item{varcov1}{
    Optional. Variance-covariance matrix matrix of the kernel with
    bandwidth \eqn{h=1}. See section on Anisotropic Smoothing.
  }
}
\details{
  This function selects an appropriate bandwidth for the nonparametric
  smoothing of mark values using \code{\link{Smooth.ppp}}.
  
  The argument \code{X} must be a marked point pattern
  with a vector or data frame of marks. All mark values must be numeric.
  
  The bandwidth is selected by least-squares cross-validation.
  Let \eqn{y_i}{y[i]} be the mark value at the \eqn{i}th data point.
  For a particular choice of smoothing bandwidth,
  let \eqn{\hat y_i}{y*[i]} be the smoothed value at the \eqn{i}th data point.
  Then the bandwidth is chosen to minimise 
  the squared error of the smoothed values
  \eqn{\sum_i (y_i - \hat y_i)^2}{sum (y[i] - y*[i])^2}.

  The result of \code{bw.smoothppp}
  is a numerical value giving the selected bandwidth \code{sigma}.
  The result also belongs to the class \code{"bw.optim"}
  allowing it to be printed and plotted. The plot shows the cross-validation
  criterion as a function of bandwidth.  
  
  The range of values for the smoothing bandwidth \code{sigma}
  is set by the arguments \code{hmin, hmax}. There is a sensible default,
  based on the nearest neighbour distances.
  
  If the optimal bandwidth is achieved at an endpoint of the
  interval \code{[hmin, hmax]}, the algorithm will issue a warning
  (unless \code{warn=FALSE}). If this occurs, then it is probably advisable
  to expand the interval by changing the arguments \code{hmin, hmax}.

  Computation time depends on the number \code{nh} of trial values
  considered, and also on the range \code{[hmin, hmax]} of values
  considered, because larger values of \code{sigma} require
  calculations involving more pairs of data points.
}
\section{Anisotropic Smoothing}{
  Anisotropic smoothing is available in \code{\link{Smooth.ppp}} using
  the argument \code{varcov} to specify the variance-covariance
  matrix of the anisotropic kernel. In order to choose the matrix
  \code{varcov}, the user can call \code{bw.smoothppp} 
  using the argument \code{varcov1} to specify a \sQuote{template}
  matrix. Scalar multiples of \code{varcov1} will be
  considered and the optimal scale factor will be determined.
  That is, \code{bw.smoothppp} will try 
  smoothing the data using \code{varcov = h^2 * varcov1} for different
  values of \code{h} ranging from \code{hmin} to \code{hmax}. 
  The result of \code{bw.smoothppp} will be the optimal value
  of the standard deviation scale factor \code{h}.
}
\value{
  A single numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  (see \code{\link[spatstat.explore]{bw.optim.object}})
  which can be plotted to show the bandwidth selection criterion
  as a function of \code{sigma}.
}
\seealso{
  \code{\link{Smooth.ppp}},
  \code{\link[spatstat.explore]{bw.optim.object}}
}
\examples{
  \testonly{op <- spatstat.options(n.bandwidth=8)}
  b <- bw.smoothppp(longleaf)
  b
  plot(b)
  \testonly{spatstat.options(op)}
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Bandwidth selection}
