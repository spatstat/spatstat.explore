\name{bw.optim.object}
\alias{bw.optim.object} %DoNotExport
\title{Class of Optimized Bandwidths}
\description{
  An object of the class \code{"bw.optim"} represents
  a tuning parameter (usually a smoothing bandwidth)
  that has been selected automatically.
  The object can be used as if it were a numerical value,
  but it can also be plotted to show the optimality criterion.
}
\details{
  An object of the class \code{"bw.optim"} represents the numerical
  value of a smoothing bandwidth, a threshold, or a similar
  tuning parameter, that has been selected by optimising 
  a criterion such as cross-validation.

  The object is a numerical value, with some attributes that retain
  information about how the value was selected.
  
  Attributes include the vector of candidate values that were examined,
  the corresponding values of the optimality criterion, the name of the
  parameter, the name of the optimality criterion, and the units in
  which the parameter is measured.

  There are methods for \code{print}, \code{plot},
  \code{summary}, \code{\link{as.data.frame}} and \code{\link{as.fv}}
  for the class \code{"bw.optim"}.

  The \code{print} method simply prints the numerical value of the
  parameter.
  The \code{summary} method prints this value, and states how
  this value was selected.

  The \code{plot} method produces a plot of the optimisation criterion
  against the candidate value of the parameter. The \code{as.data.frame}
  and \code{as.fv} methods extract this graphical information as a data
  frame or function table, respectively.
}
\seealso{
  Functions which produce objects of class \code{bw.optim} include
  \code{\link{bw.CvL}},
  \code{\link{bw.CvL.adaptive}},
  \code{\link{bw.diggle}},
  \code{\link[spatstat.linnet]{bw.lppl}},
  \code{\link{bw.pcf}},
  \code{\link{bw.ppl}},
  \code{\link{bw.relrisk}},
  \code{\link[spatstat.linnet]{bw.relrisk.lpp}},
  \code{\link{bw.smoothppp}}
  and
  \code{\link[spatstat.linnet]{bw.voronoi}}
}
\examples{
  Ns <- if(interactive()) 32 else 3
  b <- bw.ppl(redwood, srange=c(0.02, 0.07), ns=Ns)
  b
  summary(b)
  plot(b)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{attribute}
 
 
\concept{Bandwidth selection}
