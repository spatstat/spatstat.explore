\name{coxmap}
\alias{coxmap}
\alias{coxmap.ppp}
\title{
  Identify Dense and Sparse Regions in a Point Pattern
}
\description{
  Given a point pattern, classify the spatial domain into regions
  where the pattern is dense, sparse, or neither.
}
\usage{
coxmap(X, \dots)

\method{coxmap}{ppp}(X, n, alpha=0.05, \dots)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{n}{
    Maximum number of neighbours to consider. An integer greater than 1.
  }
  \item{alpha}{
    Threshold for upper and lower quantiles. A number strictly between 0 and 1.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[spatstat.geom]{as.mask}}
    to determine the pixel resolution.
  }
}
\details{
  T. F. Cox (1979) proposed a method for identifying
  dense and sparse regions in a forest stand.
  
  Each spatial location \eqn{u} is classified as
  belonging to the \dQuote{clumped} region if
  \eqn{T_n(u) < d_n}{T(u,n) < d[n]}, or 
  belonging to the \dQuote{sparse} region if
  \eqn{T_n(u) > c_n}{T(u,n) > c[n]}, and otherwise belonging
  to neither the clumped nor the sparse regions.
  Here \eqn{d_n,c_n}{d[n],c[n]} are constants depending on
  \code{n} and \code{alpha}, and
  \eqn{T_n(u)}{T(u,n)} is the statistic
  \deqn{
    T_n(u) = b_2 + b_1 \overline\lambda \pi \sum_{k=1}^n R_k(u)^2
  }{
    T(u, n) = b2 + b1 * \lambda * pi * sum[k=1,...,n] R(u, k)^2 
  }
  where \eqn{b_1,b_2}{b1, b2} are constants depending on n,
  \eqn{\overline\lambda}{\lambda} is the estimated intensity
  (number of points of \code{X} divided by area of the window of
  \code{X}), and \eqn{R_k(u)}{R(u,k)} is the distance from \eqn{u}
  to the \eqn{k}-th nearest point of \code{X}.

  This function evaluates \eqn{T_n(u)}{T(u,n)} at each pixel location
  \eqn{u} in a grid, and classifies each pixel as belonging to the
  clumped, sparse or intermediate regions. The result is a pixel image
  with categorical values giving the classification.
}
\value{
  A pixel image (object of class \code{"im"}) with factor values,
  with possible levels \code{clumped}, \code{sparse} and \code{neither}.
}
\references{
 Cox, T.F. (1979)
 A method for mapping the dense and sparse regions of a forest stand.
 \emph{Applied Statistics} \bold{28}, 14--19.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{nnclean}},
  \code{\link[spatstat.geom]{distmap}},
  \code{\link[spatstat.geom]{nnmap}},
  \code{\link[spatstat.geom]{quantess}}
}
\examples{
   plot(coxmap(redwoodfull, 5))
   points(redwoodfull, col="white", cex=0.5)
}
\keyword{spatial}
\keyword{classif}
