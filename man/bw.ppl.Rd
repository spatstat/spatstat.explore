\name{bw.ppl}
\alias{bw.ppl}
\title{
  Likelihood Cross Validation Bandwidth Selection for Kernel Density
}
\description{
  Uses likelihood cross-validation to select a smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.ppl(X, \dots, srange=NULL, ns=16, sigma=NULL, varcov1=NULL,
          weights=NULL, shortcut=TRUE, warn=TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{srange}{
    Optional numeric vector of length 2 giving the
    range of values of bandwidth to be searched.
  }
  \item{ns}{
    Optional integer giving the number of values of
    bandwidth to search.
  }
  \item{sigma}{
    Optional. Vector of values of the bandwidth to be searched.
    Overrides the values of \code{ns} and \code{srange}.
  }
  \item{varcov1}{
    Optional. Variance-covariance matrix matrix of the kernel with
    bandwidth \eqn{h=1}. See section on Anisotropic Smoothing.
  }
  \item{weights}{
    Optional. Numeric vector of weights for the points of \code{X}.
    Argument passed to \code{\link{density.ppp}}.
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link{density.ppp}}.
  }
  \item{shortcut}{
    Logical value indicating whether to speed up the calculation
    by omitting the integral term in the cross-validation criterion.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the maximum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{\sigma} is chosen to 
  maximise the point process likelihood cross-validation criterion
  \deqn{
    \mbox{LCV}(\sigma) =
    \sum_i \log\hat\lambda_{-i}(x_i) - \int_W \hat\lambda(u) \, {\rm d}u
  }{
    LCV(\sigma) = sum[i] log(\lambda[-i](x[i])) - integral[W] \lambda(u) du
  }
  where the sum is taken over all the data points \eqn{x_i}{x[i]},
  where \eqn{\hat\lambda_{-i}(x_i)}{\lambda[-i](x_i)} is the
  leave-one-out kernel-smoothing estimate of the intensity at
  \eqn{x_i}{x[i]} with smoothing bandwidth \eqn{\sigma}{\sigma},
  and \eqn{\hat\lambda(u)}{\lambda(u)} is the kernel-smoothing estimate
  of the intensity at a spatial location \eqn{u} with smoothing
  bandwidth \eqn{\sigma}{\sigma}.
  See Loader(1999, Section 5.3).

  The value of \eqn{\mbox{LCV}(\sigma)}{LCV(\sigma)} is computed
  directly, using \code{\link{density.ppp}}, 
  for \code{ns} different values of \eqn{\sigma}{\sigma}
  between \code{srange[1]} and \code{srange[2]}.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the (rescaled) mean-square error
  as a function of \code{sigma}.

  If \code{shortcut=TRUE} (the default), the computation is accelerated by
  omitting the integral term in the equation above. This is valid
  because the integral is approximately constant. 
}
\section{Anisotropic Smoothing}{
  Anisotropic kernel smoothing is available in \code{\link{density.ppp}} using
  the argument \code{varcov} to specify the variance-covariance
  matrix of the anisotropic kernel. In order to choose the matrix
  \code{varcov}, the user can call \code{bw.ppl} 
  using the argument \code{varcov1} to specify a \sQuote{template}
  matrix. Scalar multiples of \code{varcov1} will be
  considered and the optimal scale factor will be determined.
  That is, \code{bw.ppl} will try 
  smoothing the data using \code{varcov = h^2 * varcov1} for different
  values of \code{h}.
  The result of \code{bw.ppl} will be the optimal value of \code{h}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.diggle}},
  \code{\link{bw.scott}},
  \code{\link{bw.CvL}},
  \code{\link{bw.frac}}.
}
\examples{
  if(interactive()) {
    b <- bw.ppl(redwood)
    plot(b, main="Likelihood cross validation for redwoods")
    plot(density(redwood, b))
  }
  \testonly{
    b1 <- bw.ppl(redwood, srange=c(0.03, 0.07), ns=2)
    b2 <- bw.ppl(redwood, srange=c(0.03, 0.07), ns=2, shortcut=FALSE)
  }
}
\references{
  Loader, C. (1999)
  \emph{Local Regression and Likelihood}.
  Springer, New York.
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Bandwidth selection}
