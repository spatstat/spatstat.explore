\name{bw.relriskHeatppp}
\alias{bw.relriskHeatppp}
\title{
  Bandwidth Selection for Relative Risk using Diffusion
}
\description{
  Performs data-based bandwidth selection for
  the diffusion estimate of relative risk \code{\link[spatstat.explore]{relriskHeat.ppp}}
  using either likelihood cross-validation or least squares
}
\usage{
bw.relriskHeatppp(X, \dots, method = c("likelihood", "leastsquares"),
            weights = NULL, srange = NULL, ns = 16, sigma = NULL,
            leaveoneout = TRUE, verbose = TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{relriskHeat.ppp}}.
  }
  \item{method}{
    Character string specifying the cross-validation method.
    Partially matched to \code{"likelihood"} for binary likelihood
    cross-validation or \code{"leastsquares"} for least squares
    cross-validation.
  }
  \item{weights}{
    Optional numeric vector of weights associated with each point of \code{X}.
  }
  \item{srange}{
    Numeric vector of length 2 specifying a range of bandwidths to be
    considered.
  }
  \item{ns}{
    Integer. Number of candidate bandwidths to be considered.
  }
  \item{sigma}{
    Maximum smoothing bandwidth.
    A numeric value, or a pixel image, or a \code{function(x,y)}.
    Alternatively a numeric vector containing a sequence of
    candidate bandwidths.
  }
  \item{leaveoneout}{
    Logical value specifying whether intensity values at data points
    should be estimated using the leave-one-out rule.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports.
  }
}
\details{
  This algorithm selects the optimal global bandwidth for
  kernel estimation of relative risk for the dataset \code{X}
  using diffusion smoothing \code{\link[spatstat.explore]{relriskHeat}}.
  
  If \code{sigma} is a numeric value, the algorithm finds the
  optimal bandwidth \code{tau <= sigma}.

  If \code{sigma} is a pixel image or function, the algorithm
  finds the optimal fraction \code{0 < f <= 1} such that
  smoothing with \code{f * sigma} would be optimal.
}
\value{
  A numerical value giving the selected bandwidth
  (if \code{sigma} was a numeric value)
  or the selected fraction of the maximum bandwidth
  (if \code{sigma} was a pixel image or function).
  The result also belongs to the class \code{"bw.optim"} which can be
  plotted. 
}
\author{
  \adrian, \tilman and Suman Rakshit.
}
\seealso{
  \code{\link[spatstat.explore]{relriskHeat.ppp}}
}
\examples{
  ## bovine tuberculosis data
  X <- subset(btb, select=spoligotype)
  if(interactive()) {
    smax <- 40
    ns <- 16
    dimyx <- NULL
  } else {
    ## reduce data and resolution to speed up
    X <- X[c(TRUE, rep(FALSE, 7))]
    smax <- 9
    ns <- 8
    dimyx <- 32
  }
  b <- bw.relriskHeatppp(X, sigma=smax, ns=ns, dimyx=dimyx)
  b
  plot(b) 
}
\keyword{spatial}
\keyword{nonparametric}

