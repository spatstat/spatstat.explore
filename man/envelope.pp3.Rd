\name{envelope.pp3}
\alias{envelope.pp3}
\title{Simulation Envelopes of Summary Function for 3D Point Pattern}
\description{
  Computes simulation envelopes of a summary function
  for a three-dimensional point pattern.
}
\usage{
  \method{envelope}{pp3}(Y, fun=K3est, nsim=99, nrank=1, \dots, 
  funargs=list(), funYargs=funargs, simulate=NULL, verbose=TRUE, 
  transform=NULL,global=FALSE,ginterval=NULL,use.theory=NULL,
  alternative=c("two.sided", "less", "greater"),
  scale=NULL, clamp=FALSE, 
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL,
  maxnerr=nsim, rejectNA=FALSE, silent=FALSE, 
  do.pwrong=FALSE, envir.simul=NULL)
}
\arguments{
  \item{Y}{
    A three-dimensional point pattern (object of class
    \code{"pp3"}).
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a 3D point pattern. 
  }
  \item{nsim}{
    Number of simulated point patterns to be generated
    when computing the envelopes.
  }
  \item{nrank}{
    Integer. Rank of the envelope value amongst the \code{nsim} simulated
    values. A rank of 1 means that the minimum and maximum
    simulated values will be used.
  }
  \item{\dots}{
    Extra arguments passed to \code{fun}.
  }
  \item{funargs}{
    A list, containing extra arguments to be passed to \code{fun}.
  }
  \item{funYargs}{
    Optional. A list, containing extra arguments to be passed to
    \code{fun} when applied to the original data \code{Y} only.
  }
  \item{simulate}{
    Optional. Specifies how to generate the simulated point patterns.
    If \code{simulate} is an expression in the R language, then this
    expression will be evaluated \code{nsim} times,
    to obtain \code{nsim} point patterns which are taken as the
    simulated patterns from which the envelopes are computed.
    If \code{simulate} is a function, then this function will be
    repeatedly applied to the data pattern \code{Y} to obtain
    \code{nsim} simulated patterns.
    If \code{simulate} is a list of point patterns, then the entries
    in this list will be treated as the simulated patterns from which
    the envelopes are computed.
    Alternatively \code{simulate} may be an object produced by the
    \code{envelope} command: see Details.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    during the simulations.
  }
  \item{transform}{
    Optional. A transformation to be applied to the
    function values, before the envelopes are computed.
    An expression object (see Details).
  }
  \item{global}{
    Logical flag indicating whether envelopes should be pointwise
    (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
  }
  \item{ginterval}{
    Optional.
    A vector of length 2 specifying
    the interval of \eqn{r} values for the simultaneous critical
    envelopes. Only relevant if \code{global=TRUE}.
  }
  \item{use.theory}{
    Logical value indicating whether to use the theoretical value,
    computed by \code{fun}, as the reference value for simultaneous
    envelopes. Applicable only when \code{global=TRUE}.
  }
  \item{alternative}{
    Character string determining whether the envelope corresponds
    to a two-sided test (\code{side="two.sided"}, the default)
    or a one-sided test with a lower critical boundary
    (\code{side="less"}) or a one-sided test
    with an upper critical boundary (\code{side="greater"}).
  }
  \item{scale}{
    Optional. Scaling function for global envelopes.
    A function in the \R language which determines the
    relative scale of deviations, as a function of
    distance \eqn{r}, when computing the global envelopes.
    Applicable only when \code{global=TRUE}.
    Summary function values for distance \code{r}
    will be \emph{divided} by \code{scale(r)} before the
    maximum deviation is computed. The resulting global envelopes
    will have width proportional to \code{scale(r)}. 
  }
  \item{clamp}{
    Logical value indicating how to compute envelopes when
    \code{alternative="less"} or \code{alternative="greater"}.
    Deviations of the observed
    summary function from the theoretical summary function are initially
    evaluated as signed real numbers, with large positive values indicating
    consistency with the alternative hypothesis.
    If \code{clamp=FALSE} (the default), these values are not changed.
    If \code{clamp=TRUE}, any negative values are replaced by zero.
  }
  \item{savefuns}{
    Logical flag indicating whether to save all the simulated
    function values.
  }
  \item{savepatterns}{
    Logical flag indicating whether to save all the simulated
    point patterns.
  }
  \item{nsim2}{
    Number of extra simulated point patterns to be generated
    if it is necessary to use simulation to estimate the theoretical
    mean of the summary function. Only relevant when \code{global=TRUE}
    and the simulations are not based on CSR.
  }
  \item{VARIANCE}{
    Logical. If \code{TRUE}, critical envelopes will be calculated
    as sample mean plus or minus \code{nSD} times sample standard
    deviation.
  }
  \item{nSD}{
    Number of estimated standard deviations used to determine
    the critical envelopes, if \code{VARIANCE=TRUE}.
  }
  \item{Yname}{
    Character string that should be used as the name of the 
    data point pattern \code{Y} when printing or plotting the results.
  }
  \item{maxnerr}{
    Maximum number of rejected patterns.
    If \code{fun} yields a fatal error when applied to a simulated point
    pattern (for example, because the pattern is empty and \code{fun}
    requires at least one point), the pattern will be rejected
    and a new random point pattern will be generated. If this happens
    more than \code{maxnerr} times, the algorithm will give up.
  }
  \item{rejectNA}{
    Logical value specifying whether to reject a simulated pattern
    if the resulting values of \code{fun} are all equal to \code{NA},
    \code{NaN} or infinite. If \code{FALSE} (the default), then
    simulated patterns are only rejected when \code{fun} gives a
    fatal error. 
  }
  \item{silent}{
    Logical value specifying whether to print a report each time
    a simulated pattern is rejected.
  }
  \item{do.pwrong}{
    Logical. If \code{TRUE}, the algorithm will also estimate
    the true significance level of the \dQuote{wrong} test (the test that
    declares the summary function for the data to be significant
    if it lies outside the \emph{pointwise} critical boundary at any
    point). This estimate is printed when the result is printed.
  }
  \item{envir.simul}{
    Environment in which to evaluate the expression \code{simulate},
    if not the current environment.
  }
}
\value{
  A function value table (object of class \code{"fv"})
  which can be plotted directly.
  See \code{\link{envelope}} for further details.
}
\details{
  The \code{envelope} command performs simulations and
  computes envelopes of a summary statistic based on the simulations.
  The result is an object that can be plotted to display the envelopes.
  The envelopes can be used to assess the goodness-of-fit of
  a point process model to point pattern data.
  
  The \code{envelope} function is generic, with methods for
  the classes \code{"ppp"}, \code{"ppm"} and \code{"kppm"}
  described in the help file for \code{\link{envelope}}.
  This function \code{envelope.pp3} is the method for 
  three-dimensional point patterns (objects of class \code{"pp3"}).
  
  For the most basic use, if you have a 3D point pattern \code{X} and
  you want to test Complete Spatial Randomness (CSR), type
  \code{plot(envelope(X, K3est,nsim=39))} to see the three-dimensional
  \eqn{K} function for \code{X} plotted together with the envelopes of
  the three-dimensional \eqn{K} function for 39 simulations of CSR. 
  
  To create simulation envelopes, the command \code{envelope(Y, ...)} 
  first generates \code{nsim} random point patterns
  in one of the following ways. 
  \itemize{
    \item 
    If \code{simulate=NULL},
    then we generate \code{nsim} simulations of
    Complete Spatial Randomness (i.e. \code{nsim} simulated point patterns
    each being a realisation of the uniform Poisson point process)
    with the same intensity as the pattern \code{Y}.
    \item
    If \code{simulate} is supplied, then it determines how the
    simulated point patterns are generated.
    See \code{\link{envelope}} for details.
  }
  
  The summary statistic \code{fun} is applied to each of these simulated
  patterns. Typically \code{fun} is one of the functions
  \code{K3est}, \code{G3est}, \code{F3est} or \code{pcf3est}.
  It may also be a character string
  containing the name of one of these functions.
  
  For further information, see the documentation for
  \code{\link{envelope}}.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A. (1993)
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42}, 641--668.
} 
\seealso{
  \code{\link[spatstat.geom]{pp3}},
  \code{\link[spatstat.random]{rpoispp3}},
  \code{\link{K3est}},
  \code{\link{G3est}},
  \code{\link{F3est}},
  \code{\link{pcf3est}}.
}
\examples{
   X <- rpoispp3(20, box3())
 if(interactive()) {
 plot(envelope(X, nsim=39))
 }
 \testonly{
  plot(envelope(X, nsim=4))
 }
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
 
 
\concept{Three-dimensional}
\concept{Goodness-of-fit}
\concept{Test of randomness}
\concept{Envelope of simulations}
