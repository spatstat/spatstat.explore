\name{bw.CvL.adaptive}
\alias{bw.CvL.adaptive}
\title{
  Select Adaptive Bandwidth for Kernel Estimation
  Using Cronie-Van Lieshout Criterion
}
\description{
  Uses the  Cronie-Van Lieshout criterion to
  select the global smoothing bandwidth
  for adaptive kernel estimation of point process intensity.
}
\usage{
bw.CvL.adaptive(X, \dots, 
            hrange = NULL, nh = 16, h=NULL,
            bwPilot = bw.scott.iso(X),
            edge = FALSE, diggle = TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link{densityAdaptiveKernel.ppp}}.
  }
  \item{hrange}{
    Optional numeric vector of length 2 giving the
    range of values of global bandwidth \code{h} to be searched.
  }
  \item{nh}{
    Optional integer giving the number of values of
    bandwidth \code{h} to search.
  }
  \item{h}{
    Optional. Vector of values of the bandwidth to be searched.
    Overrides the values of \code{nh} and \code{hrange}.
  }
  \item{bwPilot}{
    Pilot bandwidth. A scalar value in the same units as the
    coordinates of \code{X}. The smoothing bandwidth
    for computing an initial estimate of intensity using
    \code{\link{density.ppp}}.
  }
  \item{edge}{
    Logical value indicating whether to apply edge correction.
  }
  \item{diggle}{
    Logical. If \code{TRUE}, use the Jones-Diggle improved edge correction,
    which is more accurate but slower to compute than the default
    correction.
  }
}
\details{
  This function selects an appropriate value of global bandwidth
  \code{h0} for adaptive kernel estimation of the intensity function
  for the point pattern \code{X}.

  In adaptive estimation, each point in the point pattern is
  subjected to a different amount of smoothing, controlled by
  data-dependent or spatially-varying bandwidths.
  The global bandwidth \code{h0} is a scale factor
  which is used to adjust all of the data-dependent bandwidths
  according to the Abramson (1982) square-root rule.
  
  This function considers each candidate value of bandwidth \eqn{h},
  performs the smoothing steps described above, extracts the
  adaptively-estimated intensity values
  \eqn{\hat\lambda(x_i)}{lambda(X[i])} at each data point \eqn{x_i}{X[i]},  
  and calculates the Cronie-Van Lieshout criterion
  \deqn{
    \mbox{CvL}(h) = \sum_{i=1}^n \frac 1 {\hat\lambda(x_i)}.
  }{
    CvL(h) = sum_[i=1,...n] 1/lambda(X[i]).
  }
  The value of \eqn{h} which minimises the squared difference
  \deqn{
    LP2(h) = (CvL(h) - |W|)^2
  }
  (where \code{|W|} is the area of the window of \code{X})
  is selected as the optimal global bandwidth.

  Bandwidths \code{h} are physical distance values
  expressed in the same units as the coordinates of \code{X}.
}
\value{
  A single numerical value giving the selected global bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  (see \code{\link[spatstat.explore]{bw.optim.object}})
  which can be plotted to show the bandwidth selection criterion
  as a function of \code{sigma}.
}
\references{
  Abramson, I. (1982)
  On bandwidth variation in kernel estimates --- a square root law.
  \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr

  Cronie, O and Van Lieshout, M N M (2018) A non-model-based approach to
  bandwidth selection for kernel estimators of spatial intensity functions,
  \emph{Biometrika}, \bold{105}, 455-462.

  Van Lieshout, M.N.M. (2021)
  Infill asymptotics for adaptive kernel estimators of spatial intensity.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{63} (1) 159--181.
}
\author{
  Marie-Colette Van Lieshout.
  Modified by \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{bw.optim.object}}.
  
  \code{\link{adaptive.density}},
  \code{\link{densityAdaptiveKernel.ppp}},
  \code{\link{bw.abram.ppp}},
  \code{\link{density.ppp}}.
  
  To select a \emph{fixed} smoothing bandwidth
  using the Cronie-Van Lieshout criterion, use \code{\link{bw.CvL}}.
}
\examples{
  online <- interactive()
  if(online) {
    h0 <- bw.CvL.adaptive(redwood3)
  } else {
    ## faster computation for package checker
    h0 <- bw.CvL.adaptive(redwood3, nh=8,
                          hrange=c(1/4, 4) * bw.diggle(redwood3))
  }
  plot(h0)
  plot(as.fv(h0), CvL ~ h)
  if(online) {
    Z <- densityAdaptiveKernel(redwood3, h0)
    plot(Z)
  }
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Adaptive smoothing}

\concept{Bandwidth selection}
