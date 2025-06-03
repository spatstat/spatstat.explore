\name{roc.rhohat}
\alias{roc.rhohat}
\title{
  Receiver Operating Characteristic
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern from a given intensity function.
}
\usage{
\method{roc}{rhohat}(X, \dots, high = TRUE)
}
\arguments{
  \item{X}{
    Estimate of rho function (object of class \code{"rhohat"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
}
\details{
  This command computes Receiver Operating
  Characteristic curve from an estimate of \code{\link[spatstat.explore]{rhohat}}.

  See \code{\link[spatstat.explore]{roc}} for general information about ROC curves.

  It is assumed that the rhohat function \code{X} describes the true functional
  dependency of the point process intensity on the given covariate (included in
  the \code{rhohat}-object). The ROC curve is estimated by integrating this
  function as described in \rocpapercite.
}
\value{
  Function value table (object of class \code{"fv"} and \code{"roc"})
  which can be plotted to show the ROC curve.
}
\author{
  \rocketAuthors.
}
\seealso{
\code{\link[spatstat.explore]{roc}}, \code{\link[spatstat.explore]{auc}}
}
\references{
  \rocpaper
}
\examples{
  rh <- rhohat(swedishpines, "x")
  plot(roc(rh))
}
\keyword{spatial}
