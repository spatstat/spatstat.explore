\name{bw.pcfinhom}
\alias{bw.pcfinhom}
\title{
  Cross Validated Bandwidth Selection for Inhomogeneous Pair Correlation Function
}
\description{
  Uses composite likelihood or generalized least squares 
  cross-validation to select a smoothing bandwidth
  for the kernel estimation of the inhomogeneous pair correlation function.
}
\usage{
  bw.pcfinhom(X, lambda=NULL, \dots, rmax=NULL, nr=10000, 
           cv.method=c("compLik", "leastSQ", "oracle"),
           leaveoneout=TRUE, simple=TRUE,
           fast=TRUE, srange=NULL, ns=32, use.count=TRUE,
           gtrue=NULL, 
           verbose=FALSE, warn=TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{pcfinhom}}.
  }
  \item{rmax}{
    Optional. Numeric value. Maximum value of the spatial lag distance \eqn{r} 
    for which \eqn{g(r)} should be evaluated.
  }
  \item{nr}{
    Integer. Number of subintervals for discretization of 
    [0, rmax] to use in computing numerical integrals.
  }
  \item{cv.method}{
    Choice of cross validation method: either
    \code{"compLik"}, \code{"leastSQ"} or \code{"oracle"} (partially matched).
  }
  \item{leaveoneout}{
    Logical value specifying whether to use leave-one-out estimators.
    See Details.
  }
  \item{simple}{
    Logical. Whether to use simple removal of spatial lag 
    distances. See Details. 
  }
  \item{fast}{
    Logical value indicating whether to find the optimal value
    by an optimization algorithm (\code{fast=TRUE}, the default)
    or by evaluating the objective function on an equally-spaced grid
    of bandwidth values (\code{fast=FALSE}).
  }
  \item{srange}{
    Optional. Numeric vector of length 2 giving the range of
    bandwidth values that should be searched to find the optimum
    bandwidth.
  }
  \item{ns}{
    Integer. Number of values of bandwidths at which to evaluate
    the objective function, when \code{fast=FALSE}.
  }
  \item{use.count}{
    Logical value specifying the benchmark for the calculation
    when \code{cv.method="martin"}. In this calculation,
    the sum of values \eqn{1/g(d_{i,j})}{1/g(d[i,j])} over all
    pairwise distances is compared to a benchmark.
    If \code{use.count=TRUE} (the default), the benchmark
    is simply the total number of pairs of points contributing to the sum.
    If \code{use.count=FALSE}, the theoretical expected value is
    used as the benchmark.
  }
  \item{gtrue}{
    Function in the \R language giving the true pair correlation
    function, when \code{cv.method="oracle"}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports
    during the optimization procedure.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the optimum value of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function selects an appropriate bandwidth \code{bw}
  for the kernel estimator of the pair correlation function 
  of a point process intensity computed by \code{\link{pcfinhom}}.

  With \code{cv.method="leastSQ"}, the bandwidth 
  \eqn{h} is chosen to minimise an unbiased 
  estimate of the integrated mean-square error criterion 
  \eqn{M(h)} defined in equation (4) in Guan (2007a).
  The code implements the fast algorithm of Jalilian and Waagepetersen
  (2018).
    
  With \code{cv.method="compLik"}, the bandwidth 
  \eqn{h} is chosen to maximise a likelihood 
  cross-validation criterion \eqn{CV(h)} defined in 
  equation (6) of Guan (2007b).
  
  \deqn{
    M(b) = \frac{\mbox{MSE}(\sigma)}{\lambda^2} - g(0)
  }{
    M(b) = \int_{0}^{rmax} \hat{g}^2(r;b) r dr - \sum_{u,v}
  }

  With \code{cv.method="oracle"}, the true pair correlation function
  must be provided as the argument \code{gtrue}. The bandwidth 
  \eqn{h} is chosen to minimise the integrated squared difference
  between the pcf estimate and the true pcf,
  \deqn{
    M(h) = \int_0^{\mbox{rmax}} (\hat{g}(r) - g(r))^2 dr
  }{
    M(h) = integral[0,rmax] (\hat{g}(r) - g(r))^2
  }
  
  The result is a numerical value giving the selected bandwidth.
}
\section{Definition of bandwidth}{
  The bandwidth \code{bw} returned by \code{bw.pcfinhom}
  is the standard deviation of the smoothing kernel,
  following the standard convention in \R.
  As mentioned in the documentation for
  \code{\link[stats]{density.default}} and \code{\link[spatstat.explore]{pcf.ppp}}, 
  this differs from other definitions of bandwidth that can be
  found in the literature. The scale parameter
  \code{h}, which is called the bandwidth in some literature,
  is defined differently.
  For example for the Epanechnikov kernel, \code{h} is the half-width
  of the kernel, and \code{bw=h/sqrt(5)}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{pcfinhom}}
}
\examples{
  b <- bw.pcfinhom(japanesepines)
  plot(pcfinhom(japanesepines, bw=b))
}
\references{
  \smoothpcfpaper
  
  Guan, Y. (2007a). 
  A composite likelihood cross-validation approach in selecting 
  bandwidth for the estimation of the pair correlation function. 
  \emph{Scandinavian Journal of Statistics}, 
  \bold{34}(2), 336--346.
  
  Guan, Y. (2007b). 
  A least-squares cross-validation bandwidth selection approach 
  in pair correlation function estimations. 
  \emph{Statistics & Probability Letters}, 
  \bold{77}(18), 1722--1729.

  Jalilian, A. and Waagepetersen, R. (2018)
  Fast bandwidth selection for estimation of the pair correlation
  function.
  \emph{Journal of Statistical Computation and Simulation},
  \bold{88}(10), 2001--2011.
  \url{https://www.tandfonline.com/doi/full/10.1080/00949655.2018.1428606}
}
\author{
  Rasmus Waagepetersen and Abdollah Jalilian. 
  Adapted for \pkg{spatstat} by \spatstatAuthors.
  Hacked by \adrian, \martinH and \tilman.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
