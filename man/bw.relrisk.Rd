\name{bw.relrisk}
\alias{bw.relrisk}
\alias{bw.relrisk.ppp}
\title{
  Cross Validated Bandwidth Selection for Relative Risk Estimation
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the estimation of relative risk.
}
\usage{
   bw.relrisk(X, ...)

   \method{bw.relrisk}{ppp}(X, method = "likelihood", \dots,
           nh = spatstat.options("n.bandwidth"),
           hmin=NULL, hmax=NULL, warn=TRUE)
}

\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}
    which has factor valued marks).
  }
  \item{method}{
    Character string determining the cross-validation method.
    Current options are \code{"likelihood"},
    \code{"leastsquares"} or
    \code{"weightedleastsquares"}.
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
  \item{\dots}{Additional arguments passed to \code{\link{density.ppp}}
    or to other methods for \code{bw.relrisk}.}
}
\details{
  This function selects an appropriate bandwidth for the nonparametric
  estimation of relative risk using \code{\link{relrisk}}.
  
  Consider the indicators \eqn{y_{ij}}{y[i,j]} which equal \eqn{1} when
  data point \eqn{x_i}{x[i]} belongs to type \eqn{j}, and equal \eqn{0}
  otherwise.
  For a particular value of smoothing bandwidth,
  let \eqn{\hat p_j(u)}{p*[j](u)} be the estimated
  probabilities that a point at location \eqn{u} will belong to
  type \eqn{j}. 
  Then the bandwidth is chosen to minimise either the negative likelihood,
  the squared error, or the approximately standardised squared error, of the
  indicators \eqn{y_{ij}}{y[i,j]} relative to the fitted
  values  \eqn{\hat p_j(x_i)}{p*[j](x[i])}. See Diggle (2003)
  or Baddeley et al (2015).

  The result is a numerical value giving the selected bandwidth \code{sigma}.
  The result also belongs to the class \code{"bw.optim"}
  allowing it to be printed and plotted. The plot shows the cross-validation
  criterion as a function of bandwidth.  
  
  The range of values for the smoothing bandwidth \code{sigma}
  is set by the arguments \code{hmin, hmax}. There is a sensible default,
  based on multiples of Stoyan's rule of thumb \code{\link{bw.stoyan}}.
  
  If the optimal bandwidth is achieved at an endpoint of the
  interval \code{[hmin, hmax]}, the algorithm will issue a warning
  (unless \code{warn=FALSE}). If this occurs, then it is probably advisable
  to expand the interval by changing the arguments \code{hmin, hmax}.

  Computation time depends on the number \code{nh} of trial values
  considered, and also on the range \code{[hmin, hmax]} of values
  considered, because larger values of \code{sigma} require
  calculations involving more pairs of data points.
}
\value{
  A single numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  (see \code{\link[spatstat.explore]{bw.optim.object}})
  which can be plotted to show the bandwidth selection criterion
  as a function of \code{sigma}.
}
\seealso{
  \code{\link{relrisk}},
  \code{\link{bw.stoyan}}.

  \code{\link[spatstat.explore]{bw.optim.object}}.
}
\examples{
  \testonly{op <- spatstat.options(n.bandwidth=8)}
  b <- bw.relrisk(urkiola)
  b
  plot(b)
  b <- bw.relrisk(urkiola, hmax=20)
  plot(b)
  \testonly{spatstat.options(op)}
}
\references{
  \baddrubaturnbook
  
  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
  
  Kelsall, J.E. and Diggle, P.J. (1995)
  Kernel estimation of relative risk.
  \emph{Bernoulli} \bold{1}, 3--16.
}
\author{
  \adrian and \rolf.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Bandwidth selection}
