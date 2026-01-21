\name{pcf.ppp}
\alias{pcf.ppp}
\title{Pair Correlation Function of Point Pattern}
\description{
  Estimates the pair correlation function of
  a point pattern using kernel methods.
}
\usage{
  \method{pcf}{ppp}(X, \dots, r = NULL, rmax=NULL, 
                    adaptive=FALSE,
                    kernel="epanechnikov", bw=NULL, h=NULL,
                    bw.args=list(), stoyan=0.15, adjust=1, 
                    correction=c("translate", "Ripley"),
                    divisor = c("r", "d", "a", "t"),
                    zerocor=c("weighted", "reflection", "convolution",
                             "bdrykern", "JonesFoster", "none"),
                    gref = NULL,
                    tau = 0,
                    fast = TRUE, 
                    var.approx = FALSE,
                    domain=NULL,
                    ratio=FALSE, close=NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[stats]{density.default}} or
    to \code{\link[spatstat.univar]{densityBC}} controlling the kernel smoothing
    of the pair correlation estimate.
  }
  \item{r}{
    Optional. Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
    A single numeric value.
  }
  \item{adaptive}{
    Logical value specifying whether to use adaptive kernel smoothing
    (\code{adaptive=TRUE}) or fixed-bandwidth kernel smoothing
    (\code{adaptive=FALSE}, the default).
  }
  \item{kernel}{
    Choice of smoothing kernel, passed to \code{\link[stats]{density.default}}.
  }
  \item{bw}{
    Bandwidth for smoothing kernel.
    Either a single numeric value giving the standard deviation of the
    kernel, or a character string specifying a bandwidth selection rule,
    or a function that computes the selected bandwidth.
    See Details.
  }
  \item{h}{
    Kernel halfwidth \eqn{h} (incompatible with argument \code{bw}).
    A numerical value.
    The parameter \code{h} is defined as the
    half-width of the support of the kernel, except for the Gaussian
    kernel where \code{h} is the standard deviation.
  }
  \item{bw.args}{
    Optional. List of additional arguments to be passed to \code{bw}
    when \code{bw} is a function. Alternatively, \code{bw} may be a
    function that should be applied to \code{X} to produce a list of
    additional arguments.
  }
  \item{stoyan}{
    Coefficient for Stoyan's bandwidth selection rule; see Details.
  }
  \item{adjust}{
     Numerical adjustment factor for the bandwidth.
     The bandwidth actually used is \code{adjust * bw}.
     This makes it easy to specify choices like \sQuote{half the
     selected bandwidth}.
  }
  \item{correction}{
    Edge correction. A character vector specifying the choice
    (or choices) of edge correction. See Details.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"},
    or the new alternatives \code{"a"} or \code{"t"}. See Details.
  }
  \item{zerocor}{
    String (partially matched) specifying a correction for the boundary effect
    bias at \eqn{r=0}. Possible values are
    \code{"none"}, \code{"weighted"}, \code{"convolution"},
    \code{"reflection"}, \code{"bdrykern"} and \code{"JonesFoster"}.
    See Details, or help file for \code{\link[spatstat.univar]{densityBC}}.
  }
  \item{gref}{
    Optional. A pair correlation function that will be used as the
    reference for the transformation to uniformity, when
    \code{divisor="t"}. Either a \code{function} in the \R language
    giving the pair correlation function, or a fitted model
    (object of class \code{"kppm"}, \code{"dppm"}, \code{"ppm"}
    or \code{"slrm"}) or a theoretical point process model
    (object of class \code{"zclustermodel"} or \code{"detpointprocfamily"})
    for which the pair correlation function
    can be computed.  
  }
  \item{tau}{
    Optional shrinkage coefficient. A single numeric value.
  }
  \item{fast}{
    Logical value specifying whether to compute the kernel smoothing
    using a Fast Fourier Transform algorithm (\code{fast=TRUE})
    or an exact analytic kernel sum (\code{fast=FALSE}).
  }
  \item{var.approx}{
    Logical value indicating whether to compute an analytic
    approximation to the variance of the estimated pair correlation.
  }
  \item{domain}{
    Optional. Calculations will be restricted to this subset
    of the window. See Details.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
  \item{close}{
    Advanced use only. Precomputed data. See section on Advanced Use.
  }
} 
\value{
  A function value table
  (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g(r)} for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g(r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g(r)}
    estimated by Ripley isotropic correction
  }
  \item{v}{vector of approximate values of the variance of
    the estimate of \eqn{g(r)}
  }
  as required.

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{g(r)}.

  The return value also has an attribute \code{"bw"} giving the
  smoothing bandwidth that was used, and an attribute \code{"info"}
  containing details of the algorithm parameters.
}
\details{
  The pair correlation function \eqn{g(r)} 
  is a summary of the dependence between points in a spatial point
  process. The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda^2 g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda^2 * g(r) dx dy
  }
  where \eqn{\lambda}{lambda} is the intensity of the point process.
  For a completely random (uniform Poisson) process,
  \eqn{p(r) = \lambda^2 \,{\rm d}x \, {\rm d}y}{p(r) = lambda^2 dx dy}
  so \eqn{g(r) = 1}.
  Formally, the pair correlation function of a stationary point process
  is defined by 
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}, the
  reduced second moment function (aka ``Ripley's \eqn{K} function'')
  of the point process. See \code{\link[spatstat.explore]{Kest}} for information
  about \eqn{K(r)}.

  For a stationary Poisson process, the
  pair correlation function is identically equal to 1. Values
  \eqn{g(r) < 1} suggest inhibition between points;
  values greater than 1 suggest clustering.

  This routine computes an estimate of \eqn{g(r)}
  by kernel smoothing. 

  \itemize{
    \item
    If \code{divisor="r"} (the default), then the standard
    kernel estimator (Stoyan and Stoyan, 1994, pages 284--285)
    is used. By default, the recommendations of Stoyan and Stoyan (1994)
    are followed exactly. 
    \item
    If \code{divisor="d"} then a modified estimator is used
    (Guan, 2007): the contribution from
    an interpoint distance \eqn{d_{ij}}{d[ij]} to the
    estimate of \eqn{g(r)} is divided by \eqn{d_{ij}}{d[ij]}
    instead of dividing by \eqn{r}. This usually improves the
    bias of the estimator when \eqn{r} is close to zero.
    \item
    If \code{divisor="a"} then improved method of
    \smoothpcfpapercite is used. The distances \eqn{d_{ij}}{d[ij]}
    are first converted to disc areas
    \eqn{a_{ij}=\pi d_{ij}^2}{a[ij] = pi * d[ij]^2},
    and smoothing is performed on the area scale,
    then the result is back-transformed to the original scale.
    \item
    If \code{divisor="t"} then the distances \eqn{d_{ij}}{d[ij]}
    are transformed to uniformity
    using the reference pair correlation function \code{gref}
    as described in \smoothpcfpapercite.
    \item
    If \code{divisor} is a \code{function} in the \R language, then
    it will be applied to the point pattern \code{X} and should return
    one of the strings \code{"r"}, \code{"d"}, \code{"a"} or \code{"t"}
    listed above. This option makes it possible to specify a rule
    which decides which estimator to use, based on the data.
  }

  There is also a choice of spatial edge corrections
  (which are needed to avoid bias due to edge effects
  associated with the boundary of the spatial window):

  \itemize{
    \item
    If \code{correction="translate"} or \code{correction="translation"}
    then the translation correction
    is used. For \code{divisor="r"} the translation-corrected estimate
    is given in equation (15.15), page 284 of Stoyan and Stoyan (1994).
    \item
    If \code{correction="Ripley"} or \code{correction="isotropic"}
    then Ripley's isotropic edge correction
    is used. For \code{divisor="r"} the isotropic-corrected estimate
    is given in equation (15.18), page 285 of Stoyan and Stoyan (1994). 
    \item
    If \code{correction="none"} then no edge correction is used,
    that is, an uncorrected estimate is computed. 
  }
  Multiple corrections can be selected. The default is
  \code{correction=c("translate", "Ripley")}.
  
  Alternatively \code{correction="all"} selects all options;
  \code{correction="best"} selects the option which has the best
  statistical performance; \code{correction="good"} selects the option
  which is the best compromise between statistical performance and speed
  of computation.

  Argument \code{zerocor} determines the correction
  to the one-dimensional kernel-smoothed estimate
  on the real number line, to correct bias close to the boundary \eqn{r=0}.
  The argument \code{zerocor} is passed to 
  \code{\link[spatstat.univar]{densityBC}}.
  Options include:
  \itemize{
    \item
    \code{zerocor="none"}: no correction.
    \item
    \code{zerocor="convolution"}: the convolution, uniform or
     renormalization kernel.
    \item
    \code{zerocor="weighted"}: the cut-and-normalization method.
    \item
    \code{zerocor="reflection"}:
    the reflection method.
    \item
    \code{zerocor="bdrykern"}: the linear boundary kernel.
    \item
    \code{zerocor="JonesFoster"}:
      the Jones-Foster modification of the linear boundary kernel.
  }

  The choice of smoothing kernel is controlled by the 
  argument \code{kernel} which is passed
  to \code{\link[stats]{density.default}}.
  The default is the Epanechnikov kernel, recommended by
  Stoyan and Stoyan (1994, page 285).

  The bandwidth of the smoothing kernel can be controlled by the
  argument \code{bw}. Bandwidth is defined as the standard deviation
  of the kernel; see the documentation for \code{\link[stats]{density.default}}.
  For the Epanechnikov kernel with half-width \code{h},
  the argument \code{bw} is equivalent to \eqn{h/\sqrt{5}}{h/sqrt(5)}.

  Stoyan and Stoyan (1994, page 285) recommend using the Epanechnikov
  kernel with support \eqn{[-h,h]} chosen by the rule of thumn
  \eqn{h = c/\sqrt{\lambda}}{h = c/sqrt(lambda)},
  where \eqn{\lambda}{lambda} is the (estimated) intensity of the
  point process, and \eqn{c} is a constant in the range from 0.1 to 0.2.
  See equation (15.16).
  If \code{bw} is missing or \code{NULL},
  then this rule of thumb will be applied.
  The argument \code{stoyan} determines the value of \eqn{c}.
  The smoothing bandwidth that was used in the calculation is returned
  as an attribute of the final result.

  The argument \code{bw} can be
  \itemize{
    \item
    missing or null. In this case, the default value for \code{bw}
    is \code{"stoyan"} when \code{adaptive=FALSE}
    and \code{"bw.abram"} when \code{adaptive=TRUE}.
    \item
    a single numeric value giving the bandwidth.
    \item
    a character string specifying a bandwidth selection rule.
    String names of rules applicable when \code{adaptive=FALSE} 
    include \code{"stoyan"}, \code{"fiksel"}
    and any rules recognised by \code{\link[stats]{density.default}}.
    String names applicable when \code{adaptive=TRUE} include
    \code{"bw.abram"} and \code{"bw.pow"}.
    \item
    a function that computes the selected bandwidth.
    \itemize{
      \item 
      If \code{adaptive=FALSE}, the function \code{bw} will be applied to the
      point pattern \code{X} to determine the bandwidth. Examples include
      \code{\link{bw.pcf}} and \code{\link{bw.stoyan}}.
      The function \code{bw} should accept the point pattern \code{X}
      as its first argument. Additional arguments to \code{bw} may be
      specified in the list \code{bw.args}. If \code{bw} recognises any
      of the arguments
      \code{kernel}, \code{correction}, \code{divisor}, \code{zerocor}
      and \code{adaptive}, then these arguments will be passed to
      \code{bw} as well. 
      The function \code{bw} should return a single
      numeric value.
      \item
      If \code{adaptive=TRUE}, the function \code{bw} will be applied to
      the vector of pairwise distances between data points (or the
      transformed distances if \code{divisor="a"} or \code{divisor="t"}).
      Examples include \code{\link[spatstat.univar]{bw.abram.default}}
      and \code{\link[spatstat.univar]{bw.pow}}.
      The function \code{bw} should accept the vector of pairwise
      distances as its first argument. Additional arguments to \code{bw} may be
      specified in the list \code{bw.args}. 
     }
  }
  Note that if \code{bw.args} is a function, it will be applied to
  the point pattern \code{X} to determine the list of arguments
  (whether \code{adaptive} is \code{TRUE} or \code{FALSE}).

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{g(r)} should be evaluated.
  There is a sensible default.
  If it is specified, \code{r} must be a vector of increasing numbers
  starting from \code{r[1] = 0}, 
  and \code{max(r)} must not exceed half the diameter of 
  the window.

  If the argument \code{domain} is given, estimation will be restricted
  to this region. That is, the estimate of 
  \eqn{g(r)} will be based on pairs of points in which the first point lies
  inside \code{domain} and the second point is unrestricted.
  The argument \code{domain}
  should be a window (object of class \code{"owin"}) or something acceptable to
  \code{\link[spatstat.geom]{as.owin}}. It must be a subset of the
  window of the point pattern \code{X}.

  To compute a confidence band for the true value of the
  pair correlation function, use \code{\link[spatstat.explore]{lohboot}}.

  If \code{var.approx = TRUE}, the variance of the
  estimate of the pair correlation will also be calculated using
  an analytic approximation (Illian et al, 2008, page 234)
  which is valid for stationary point processes which are not
  too clustered. This calculation is not yet implemented when
  the argument \code{domain} is given.

  If \code{fast=TRUE}, the calculation uses the Fast Fourier Transform
  to the maximum extent possible for the chosen boundary correction.
  If \code{fast=FALSE} (the default), the entire calculation uses
  analytic formulas written in \code{C} code.
}
\section{Advanced Use}{
  To perform the same computation using several different bandwidths \code{bw},
  it is efficient to use the argument \code{close}.
  This should be the result of \code{\link[spatstat.geom]{closepairs}(X, rmax)}
  for a suitably large value of \code{rmax}, namely
  \code{rmax >= max(r) + 3 * bw}.
}
\references{
  \smoothpcfpaper
  
  Guan, Y. (2007)
  A least-squares cross-validation bandwidth selection approach in pair
  correlation function estimation.
  \emph{Statistics and Probability Letters} \bold{77} (18) 1722--1729.

  Illian, J., Penttinen, A., Stoyan, H. and Stoyan, D. (2008)
  \emph{Statistical Analysis and Modelling of Spatial Point Patterns.}
  Wiley.

  Stoyan, D. and Stoyan, H. (1994)
  \emph{Fractals, random shapes and point fields:
  methods of geometrical statistics.}
  John Wiley and Sons.
}
\seealso{
  \code{\link[spatstat.univar]{densityBC}},
  
  \code{\link[spatstat.explore]{Kest}},
  \code{\link[spatstat.explore]{pcf}},
  \code{\link[stats]{density.default}},
  \code{\link[spatstat.explore]{bw.stoyan}},
  \code{\link[spatstat.explore]{bw.pcf}},
  \code{\link[spatstat.explore]{lohboot}}.
}
\examples{
  pr <- pcf(redwood, divisor="r")
  plot(pr, main="pair correlation function for redwoods")

  # compare estimates
  pd <- pcf(redwood, divisor="d")
  pa <- pcf(redwood, divisor="a")

  plot(pr, cbind(iso, theo) ~ r, col=c("red", "black"),
       ylim.covers=0, main="Estimates of PCF",
       lwd=c(2,1), lty=c(1,3), legend=FALSE)
  plot(pd, iso ~ r, col="blue", lwd=2, add=TRUE)
  plot(pa, iso ~ r, col="green", lwd=2, add=TRUE)
  legend("center", col=c("red", "blue", "green"), lty=1, lwd=2,
         legend=c("divisor=r","divisor=d", "divisor=a"))
}
\author{
  \spatstatAuthorsComma, \martinH  and \tilman.
}
\keyword{spatial}
\keyword{nonparametric}
