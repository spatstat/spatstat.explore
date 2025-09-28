\name{nncount}
\alias{nncount}
\title{
  Nearest Neighbour Cross-Type Occurrence Count Function
}
\description{
  Compute the estimated nearest neighbour cross-type occurrence count function
  for a multitype point pattern dataset.
}
\usage{
nncount(X, i = 1, j = 2, \dots, kmax = 20, ratio = TRUE, cumulative = TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern
    (object of class \code{"ppp"} with marks that are a factor).
  }
  \item{i}{The type (mark value)
    of the points in \code{X} which will serve as the reference points
    from which neighbours are counted.
    A character string (or something that will be converted to a
    character string) or an integer.
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} which will be counted as neighbours of interest.
    A character string (or something that will be
    converted to a character string), or an integer.
    Defaults to the second level of \code{marks(X)}.
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
  The nearest-neighbour cross-type occurrence function
  counts the frequency with which points of types \code{i} and \code{j}
  are found close together, as defined by the nearest-neighbour 
  relationship of order \code{k} for \code{k=1,2,...,kmax}.

  The algorithm visits each point of \code{X} that belongs to
  type \code{i}, and finds its nearest neighbour, second-nearest
  neighbour, and so on, until the nearest neighbour of order
  \code{kmax}.
  The algorithm counts the number of occurrences of 
  type \code{j}.
  
  The result is a function that can be plotted to show the
  observed and expected frequencies of these occurrences as a
  function of the nearest-neighbour order \code{k}.

  The function has two versions: cumulative
  and non-cumulative.
  The non-cumulative version is a function
  \eqn{P_{ij}(k)}{P[i,j](k)} which gives, for each integer \eqn{k},
  the estimated probability that a typical point in the point pattern
  of type \code{i} will have a \eqn{k}-th nearest neighbour of type \code{j}.
  For example \eqn{P_{A,B}(5)=0.6}{P[A,B](k)} would mean that,
  for a typical point of type \code{A}, there is a 60 percent chance
  that the fifth-nearest neighbour will be a point of type \code{B}.
  When \code{cumulative=FALSE}, the algorithm computes an estimate
  of \eqn{P_{ij}(k)}{P[i,j](k)} for each integer \eqn{k} up to \code{kmax}, and
  returns the result as a function.

  The cumulative version is
  \eqn{Q_{ij}(k) = \sum_{m=1}^k P_{ij}(m)/k}{Q[i,j](k) = sum[m=1,...,k] P[i,j](m)/k}.
  This means that \eqn{Q_{ij}(k)}{Q[i,j](k)} is the average proportion
  of points of type \code{j}
  amongst the nearest \eqn{k} neighbours of a point of type \code{i}.
  For example \eqn{Q_{A,B}(5)=0.6}{Q[A,B](5)=0.6} would mean that, on average,
  three out of the first five nearest neighbours of a point of type \code{A}
  will have type \code{B}.
  This function is computed when \code{cumulative=TRUE} (the default).

  In either case, the result is a function value table (class \code{"fv"})
  with two columns of values, one giving the estimate of
  \eqn{P_{i,j}(k)}{P[i,j](k)} or \eqn{Q_{i,j}(k)}{Q[i,j](k)},
  and the other giving the corresponding
  average value for the pattern. Estimated values greater than the average
  suggest that points of type \code{j} are clustered closer to
  points of type \code{i} than would be expected by chance.

  The arguments \code{i} and \code{j} should specify levels of the factor
  \code{marks(X)}. They can be given either as character strings
  (which will be matched to the levels of the factor)
  or integers (which will be interpreted as indices for the
  vector of levels of the factor, so that \code{1} represents the first level.)
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
\code{\link[spatstat.explore]{nnequal}}
}
\examples{
  plot(nncount(amacrine, "off", "on"))
  plot(envelope(amacrine, nncount, nsim=19))
  plot(envelope(amacrine, nncount, cumulative=FALSE, nsim=19))
}
\keyword{spatial}
\keyword{nonparametric}
