\name{clarkevans.test}
\alias{clarkevans.test}
\title{Clark and Evans Test}
\description{
  Performs the Clark-Evans test of aggregation
  for a spatial point pattern.
}
\usage{
clarkevans.test(X, ...,
               correction,
               clipregion=NULL,
               alternative=c("two.sided", "less", "greater",
                             "clustered", "regular"),
               method=c("asymptotic", "MonteCarlo"),
               nsim=999)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Ignored.}
  \item{correction}{
    Character string.
    The type of edge correction to be applied.
    See \code{\link{clarkevans}} and Details below.
  }
  \item{clipregion}{
    Clipping region for the guard area correction.
    A window (object of class \code{"owin"}).
    See \code{\link{clarkevans}}
  }
  \item{alternative}{
    String indicating the type of alternative for the
    hypothesis test. Partially matched.
  }
  \item{method}{
    Character string (partially matched) specifying how to calculate
    the \eqn{p}-value of the test. See Details.
  }
  \item{nsim}{
    Number of Monte Carlo simulations to perform, if a Monte Carlo
    \eqn{p}-value is required.
  }
}
\details{
  This command uses the Clark and Evans (1954) aggregation index \eqn{R}
  as the basis for a crude test of clustering or ordering of a point pattern.
  
  The Clark-Evans \emph{aggregation index} \eqn{R} is computed by the separate
  function \code{\link{clarkevans}}.

  This command \code{clarkevans.text}
  performs a hypothesis test of clustering or ordering of
  the point pattern \code{X} based on the Clark-Evans index \eqn{R}.
  The null hypothesis is Complete
  Spatial Randomness, i.e.\ a uniform Poisson process. The alternative
  hypothesis is specified by the argument \code{alternative}:
  \itemize{
    \item \code{alternative="less"} or \code{alternative="clustered"}:
    the alternative hypothesis
    is that \eqn{R < 1} corresponding to a clustered point pattern;
    \item \code{alternative="greater"} or \code{alternative="regular"}:
    the alternative hypothesis
    is that \eqn{R > 1} corresponding to a regular or ordered point pattern;
    \item \code{alternative="two.sided"}:
    the alternative hypothesis is that \eqn{R \neq 1}{R != 1}
    corresponding to a clustered or regular pattern.
  }

  The Clark-Evans index \eqn{R} is first computed for the point pattern
  dataset \code{X} using the edge correction determined by
  the arguments \code{correction} and \code{clipregion}. These arguments
  are documented in the help file for \code{\link{clarkevans}}.

  If \code{method="asymptotic"} (the default),
  the \eqn{p}-value for the test is computed by standardising
  \eqn{R} as proposed by Clark and Evans (1954) and referring the
  standardised statistic to the standard Normal distribution. 
  For this asymptotic test, the default edge correction is
  \code{correction="Donnelly"} if the window of \code{X} is a rectangle,
  and \code{correction="cdf"} otherwise. It is strongly recommended
  to avoid using \code{correction="none"} which would lead to a severely
  biased test.

  If \code{method="MonteCarlo"}, the \eqn{p}-value for the test is computed
  by comparing the observed value of \eqn{R} to the
  results obtained from \code{nsim} simulated realisations of
  Complete Spatial Randomness conditional on the
  observed number of points. This test is theoretically exact
  for any choice of edge correction, but may have lower power
  than the asymptotic test.
  For this Monte Carlo test, the default edge correction
  is \code{correction="none"} for computational efficiency.
}
\value{
  An object of class \code{"htest"} representing the result of the test.
}
\references{
  Clark, P.J. and Evans, F.C. (1954)
  Distance to nearest neighbour as a measure of spatial
  relationships in populations. \emph{Ecology} \bold{35},
  445--453.
  
  Donnelly, K. (1978) Simulations to determine the variance
  and edge-effect of total nearest neighbour distance.
  In \emph{Simulation methods in archaeology},
  Cambridge University Press, pp 91--95.
}
\author{
  \adrian
}
\seealso{
  \code{\link{clarkevans}},
  \code{\link{hopskel.test}}
}
\examples{
  # Redwood data - clustered
  clarkevans.test(redwood)
  clarkevans.test(redwood, alternative="clustered")
  clarkevans.test(redwood, correction="cdf", method="MonteCarlo", nsim=39)
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{htest}
\concept{Test of randomness}
