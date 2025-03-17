\name{bw.stoyan}
\alias{bw.stoyan}
\title{
  Stoyan's Rule of Thumb for Bandwidth for Estimating Pair Correlation
}
\description{
  Computes a rough estimate of the appropriate bandwidth
  for kernel smoothing estimators of the pair correlation function.
}
\usage{
bw.stoyan(X, co=0.15, extrapolate=FALSE, \dots)
}

\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{co}{
    Coefficient appearing in the rule of thumb. See Details.
  }
  \item{extrapolate}{
    Logical value specifying whether to use the extrapolated version
    of the rule. See Details.
  }
  \item{\dots}{Ignored.}
}
\details{
  Estimation of the pair correlation function (and similar quantities)
  by smoothing methods requires a choice of the smoothing bandwidth.
  Stoyan and Stoyan (1995, equation (15.16), page 285) proposed a
  rule of thumb for choosing the smoothing bandwidth.
  
  For the Epanechnikov kernel, the rule of thumb is to set
  the kernel's half-width \eqn{h} to
  \eqn{0.15/\sqrt{\lambda}}{0.15/sqrt(\lambda)} where
  \eqn{\lambda}{\lambda} is the estimated intensity of the point pattern,
  typically computed as the number of points of \code{X} divided by the
  area of the window containing \code{X}.

  For a general kernel, the corresponding rule is to set the
  standard deviation of the kernel to
  \eqn{\sigma = 0.15/\sqrt{5\lambda}}{\sigma = 0.15/sqrt(5 * \lambda)}.
  
  The coefficient \eqn{0.15} can be tweaked using the
  argument \code{co}.

  To ensure the bandwidth is finite, an empty point pattern is treated
  as if it contained 1 point.

  The original version of Stoyan's rule, stated above, was developed
  by experience with patterns of 30 to 100 points. For patterns with
  larger numbers of points, the bandwidth should be smaller:
  the theoretically optimal bandwidth
  decreases in proportion to \eqn{n^{-1/5}}{n^(-1/5)} where
  \eqn{n} is the number of points in the pattern.
  In the \sQuote{extrapolated} version of Stoyan's rule proposed
  by \smoothpcfpapercite, the value \eqn{\sigma} calculated above
  is multiplied by \eqn{(100/n)^{1/5}}{(100/n)^(1/5)}.
  The extrapolated rule is applied if \code{extrapolate=TRUE}.
}
\value{
  A finite positive numerical value giving the selected bandwidth (the standard
  deviation of the smoothing kernel).
}
\seealso{
  \code{\link{pcf}},
  \code{\link{bw.relrisk}}
}
\examples{
  bw.stoyan(shapley)
  bw.stoyan(shapley, extrapolate=TRUE)
}
\references{
  \smoothpcfpaper
  
  Stoyan, D. and Stoyan, H. (1995)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\author{\adrian, \rolf, \tilman, \martinH and \yamei.}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Bandwidth selection}
