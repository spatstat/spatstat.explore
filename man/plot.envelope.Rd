\name{plot.envelope}
\alias{plot.envelope}
\title{Plot a Simulation Envelope}
\description{
  Plot method for the class \code{"envelope"}.
}
\usage{
 \method{plot}{envelope}(x, \dots, main)
}
\arguments{
  \item{x}{
    An object of class \code{"envelope"}, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{main}{Main title for plot.}
  \item{\dots}{
    Extra arguments passed to \code{\link{plot.fv}}.
  }
}
\value{
  Either \code{NULL}, or a data frame giving the meaning of the
  different line types and colours.
}
\details{
  This is the \code{plot} method for the class \code{"envelope"}
  of simulation envelopes. Objects of this class are
  created by the command \code{\link{envelope}}.

  This plot method is currently identical to \code{\link{plot.fv}}.

  Its default behaviour is to shade the region
  between the upper and lower envelopes in a light grey colour.
  To suppress the shading and plot the upper and lower envelopes
  as curves, set \code{shade=NULL}.
  To change the colour of the shading, use the argument \code{shadecol}
  which is passed to \code{\link{plot.fv}}. 
  
  See \code{\link{plot.fv}} for further information on how to
  control the plot.
}
\examples{
   E <- envelope(cells, Kest, nsim=19)
   plot(E)
   plot(E, sqrt(./pi) ~ r)
}
\seealso{
  \code{\link{envelope}},
  \code{\link{plot.fv}}
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{hplot}
 
\concept{Goodness-of-fit}
\concept{Test of randomness}
\concept{Envelope of simulations}
