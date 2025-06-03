\name{youden}
\alias{youden}
\title{
  Youden Statistic
}
\description{
  Calculate the Youden statistic for an ROC curve.
}
\usage{
youden(X, sign = c("positive", "absolute", "negative"))
}
\arguments{
  \item{X}{
    ROC curve (object of class \code{"roc"} produced by
    \code{\link[spatstat.explore]{roc}}).
  }
  \item{sign}{
    Character string indicating which version of the statistic to calculate.
  }
}
\details{
  For a receiver operating characteristic (ROC) curve, the
  Youden statistic is the maximum vertical deviation between the curve 
  and the diagonal line \eqn{y=x}.

  Suppose \eqn{R(p)} denotes the ROC curve as a function of the
  horizontal coordinate \eqn{p}.
  If \code{sign="positive"} (the default), deviation is defined as the
  positive part \eqn{\max(0, R(p)-p)}{max(0, R(p)-p)}.
  If \code{sign="absolute"}, deviation is defined as the
  absolute value \eqn{| R(p) - p|}{abs(R(p) - p)}.
  If \code{sign="negative"}, deviation is defined as the
  negative part \eqn{\max(0, p - R(p))}{max(0, p - R(p))}.
  The maximum deviation over all values of \eqn{p} is determined.
  The result is always nonnegative.
}
\value{
  Numeric value or vector, containing nonnegative values.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{roc}},
  \code{\link[spatstat.explore]{auc}}
}
\references{
  Youden, W.J. (1950)
  Index for rating diagnostic tests.
  \emph{Cancer} \bold{3}, 32--35.
}
\examples{
  a <- with(split(mucosa), roc(ECL, "y", baseline=other, high=FALSE))
  youden(a)
}
\keyword{nonparametric}
