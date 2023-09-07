\name{compileCDF}
\alias{compileCDF}
\title{
  Generic Calculation of Cumulative Distribution Function of Distances
}
\description{
  A low-level function which 
  calculates the estimated cumulative distribution function
  of a distance variable.
}
\usage{
compileCDF(D, B, r, \dots, han.denom=NULL, check=TRUE)
}
\arguments{
  \item{D}{
    A vector giving the distances from each data point to the target.
  }
  \item{B}{
    A vector giving the distances from each data point to the window boundary,
    or censoring distances.
  }
  \item{r}{
    An equally spaced, finely spaced sequence of distance values
    at which the CDF should be estimated.
  }
  \item{\dots}{
    Ignored.
  }
  \item{han.denom}{
    Denominator for the Hanisch-Chiu-Stoyan estimator.
    A single number, or a numeric vector with the same length
    as \code{r}. 
  }
  \item{check}{
    Logical value specifying whether to check validity of the data,
    for example, that the vectors \code{D} and \code{B} have the same
    length, and contain non-negative numbers.
  }
}
\details{
  This low-level function calculates estimates of
  the cumulative distribution function
  \deqn{F(r) = P(D \le r)}{F(r) = P(D <= r)}
  of a distance variable \eqn{D}, given a vector of observed values of
  \eqn{D} and other information.
  Examples of this concept include the empty space distance function
  computed by \code{\link{Fest}} and the nearest-neighbour distance
  distribution function \code{\link{Gest}}.

  This function \code{compileCDF}
  and its siblings \code{\link{compileK}} and \code{\link{compilepcf}}
  are useful for code development and for teaching,
  because they perform a common task, and do the housekeeping required to
  make an object of class \code{"fv"} that represents the estimated
  function. However, they are not very efficient. 

  The argument \code{D} should be a numeric vector of shortest distances
  measured from each \sQuote{query} point to the \sQuote{target} set.
  The argument \code{B} should be a numeric vector of shortest distances
  measured from each \sQuote{query} point to the boundary of the window
  of observation.
  All entries of \code{D} and \code{B} should be non-negative.

  \code{compileCDF} calculates estimates of the cumulative distribution
  function \eqn{F(r)} using the border method (reduced sample
  estimator), the Kaplan-Meier estimator and, if \code{han.denom} is
  given, the Hanisch-Chiu-Stoyan estimator.
  See Chapter 8 of Baddeley, Rubak and Turner (2015).

  The result is an object of class \code{"fv"} representing the
  estimated function.
  Additional columns (such as a column giving the theoretical
  value) must be added by the user, with the aid of
  \code{\link{bind.fv}}.
}
\value{
  An object of class \code{"fv"} representing the estimated function.
}
\author{
  \adrian
}
\seealso{
  \code{\link{compileK}}.

  \code{\link{bind.fv}} to add more columns.
}
\references{
  \baddrubaturnbook
}
\examples{
  ## Equivalent to Gest(japanesepines)
  X <- japanesepines
  D <- nndist(X)
  B <- bdist.points(X)
  r <- seq(0, 0.25, by=0.01)
  H <- eroded.areas(Window(X), r)
  G <- compileCDF(D=D, B=B, r=r, han.denom=H)
  G <- rebadge.fv(G, new.fname="G", new.ylab=quote(G(r)))
  plot(G)
}
\keyword{spatial}
\keyword{nonparametric}
