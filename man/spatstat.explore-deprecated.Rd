\name{spatstat.explore-deprecated}
\alias{spatstat.explore-deprecated} %DoNotExport
\alias{evalCovar}
\alias{which.max.im}
\title{Deprecated spatstat.explore functions}
\description{
  Deprecated spatstat.explore functions.
}
\usage{
evalCovar(model, covariate, \dots)
which.max.im(x)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat.explore} package.
  

  \code{which.max.im(x)} is replaced by
  \code{\link[spatstat.geom]{im.apply}(x, which.max)}.
  
  The internal function \code{evalCovar}
  is replaced by the internal function \code{spatialCovariateEvidence}.
  
}
\value{
\code{which.max.im} returns an integer.
}
\keyword{internal}
