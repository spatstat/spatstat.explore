\name{nnequal}
\alias{nnequal}
\title{
  Nearest Neighbour Type-Equality Function
}
\description{
  Compute the estimated nearest neighbour type-equality function
  for a multitype point pattern dataset.
}
\usage{
nnequal(X, \dots, kmax = 20, ratio = TRUE, cumulative = TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern
    (object of class \code{"ppp"} with marks that are a factor).
  }
  \item{kmax}{
    Maximum number of neighbours to be considered for each point.
  }
  \item{cumulative}{
    Logical value indicating whether to compute the cumulative estimate.
    See Details.
  }
  \item{ratio}{
    Logical value indicating whether to save ratio information
    for future use in pooling data.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The nearest-neighbour type-equality function has two versions: cumulative
  and non-cumulative.

  The non-cumulative version is a function
  \eqn{P(k)} which gives, for each integer \eqn{k},
  the probability that a typical point in the point pattern
  will have the same type as its \eqn{k}-th neighbour.
  For example \eqn{P(5)=0.6} would mean that there is a
  60 percent chance that the fifth-nearest neighbour of a point
  has the same type as the original point.
  When \code{cumulative=FALSE}, the algorithm computes an estimate
  of \eqn{P(k)} for each integer \eqn{k} up to \code{kmax}, and
  returns the result as a function.
  This is equivalent to
  computing the nearest-neighbour correlation \code{\link[spatstat.explore]{nncorr}(X, k)}
  for each integer \code{k} from \code{1} to \code{kmax}.

  The cumulative version is
  \eqn{Q(k) = \sum_{m=1}^k P(m)/k}{Q(k) = sum[m=1,...,k] P(m)/k}.
  This means that \eqn{Q(k)} is the average proportion,
  amongst the nearest \eqn{k} neighbours of a point,
  which share the same type as the original point.
  For example \eqn{Q(5)=0.6} would mean that, on average,
  three out of the first five nearest neighbours of a point
  have the same type as the original point.
  This function is computed when \code{cumulative=TRUE} (the default).

  In either case, the result is a function value table (class \code{"fv"})
  with two columns of values, one giving the estimate of
  \eqn{P(k)} or \eqn{Q(k)}, and the other giving the corresponding
  average value for the pattern. Estimated values greater than the average
  suggest that points of the same type are clustered together.
}
\value{
  A function value table (object of class \code{"fv"}) which can be
  printed and plotted. The function has argument \code{k} (the order of
  neighbour) and columns of values labelled \code{bord} (for the
  border-correction estimate) and \code{theo} (for the theoretical value
  expected under random labelling).
}
\references{
  Diez-\Martin, F., Cobo-\Sanchez, L., Baddeley, A., Uribelarrea, D.,
  Mabulla, A., Baquedano, E. and \Dominguez-Rodrigo, M. (2021)
  Tracing the spatial imprint of Oldowan technological behaviors: A view
  from DS (Bed I, Olduvai Gorge, Tanzania).
  \emph{PLOS ONE, Public Library of Science}, \bold{16}, 1--47.
  \code{DOI: 10.1371/journal.pone.0254603}
}
\author{
  \adrian and \Lucia Cobo-\Sanchez.
}
\seealso{
  \code{\link[spatstat.explore]{nncount}}
}
\examples{
  plot(nnequal(lansing), ylim=c(0, 0.5))
  plot(envelope(amacrine, nnequal, nsim=19))
  plot(envelope(amacrine, nnequal, cumulative=FALSE, nsim=19))
}
\keyword{spatial}
\keyword{nonparametric}
