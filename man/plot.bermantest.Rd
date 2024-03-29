\name{plot.bermantest}
\alias{plot.bermantest}
\title{Plot Result of Berman Test}
\description{
  Plot the result of Berman's test of goodness-of-fit
}
\usage{
\method{plot}{bermantest}(x, ...,
                   lwd=par("lwd"), col=par("col"), lty=par("lty"),
                   lwd0=lwd, col0=2, lty0=2)
}
\arguments{
  \item{x}{
    Object to be plotted. An object of class \code{"bermantest"}
    produced by \code{\link{berman.test}}.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting function
    \code{\link{plot.ecdf}}.
  }
  \item{col,lwd,lty}{
    The width, colour and type of lines used to plot the
    empirical distribution curve.
  }
  \item{col0,lwd0,lty0}{
    The width, colour and type of lines used to plot the
    predicted (null) distribution curve.
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for the class \code{"bermantest"}.
  An object of this class represents the outcome of Berman's test
  of goodness-of-fit of a spatial Poisson point process model,
  computed by \code{\link{berman.test}}.

  For the \emph{Z1} test (i.e. if \code{x} was computed using
  \code{berman.test( ,which="Z1")}), 
  the plot displays the two cumulative distribution functions
  that are compared by the test: namely the empirical cumulative distribution
  function of the covariate at the data points, \eqn{\hat F}{Fhat},
  and the predicted
  cumulative distribution function of the covariate under the model,
  \eqn{F_0}{F0}, both plotted against the value of the covariate.
  Two vertical lines show the mean values of these two distributions.
  If the model is correct, the two curves should be close; the test is
  based on comparing the two vertical lines.

  For the \emph{Z2} test (i.e. if \code{x} was computed using
  \code{berman.test( ,which="Z2")}), the plot displays the empirical
  cumulative distribution function of the values
  \eqn{U_i = F_0(Y_i)}{U[i] = F0(Y[i])} where \eqn{Y_i}{Y[i]} is the
  value of the covariate at the \eqn{i}-th data point. The diagonal line
  with equation \eqn{y=x} is also shown. Two vertical lines show the
  mean of the values \eqn{U_i}{U[i]} and the value \eqn{1/2}. If the
  model is correct, the two curves should be close. The test is based on
  comparing the two vertical lines. 
}
\seealso{
  \code{\link{berman.test}}
}
\examples{
   plot(berman.test(cells, "x"))

   if(require("spatstat.model")) {
     # synthetic data: nonuniform Poisson process
     X <- rpoispp(function(x,y) { 100 * exp(-x) }, win=square(1))

     # fit uniform Poisson process
     fit0 <- ppm(X ~1)

     # test covariate = x coordinate
     xcoord <- function(x,y) { x }

     # test wrong model
     k <- berman.test(fit0, xcoord, "Z1")
   
     # plot result of test
     plot(k, col="red", col0="green")

     # Z2 test
     k2 <- berman.test(fit0, xcoord, "Z2")
     plot(k2, col="red", col0="green")
   }
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
\concept{Goodness-of-fit}
