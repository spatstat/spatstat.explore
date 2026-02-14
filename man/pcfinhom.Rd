\name{pcfinhom}
\alias{pcfinhom}
\title{
  Inhomogeneous Pair Correlation Function
}
\description{
  Estimates the inhomogeneous pair correlation function of
  a point pattern using kernel methods.
}
\usage{
pcfinhom(X, lambda = NULL, \dots,
            r = NULL, rmax=NULL, adaptive = FALSE,
            kernel = "epanechnikov", bw = NULL, h = NULL,
            bw.args = list(), stoyan = 0.15, adjust = 1,
            correction = c("translate", "Ripley"),
            divisor=c("a", "r", "d", "t"),
            zerocor=c("convolution", "reflection", "bdrykern",
                      "JonesFoster", "weighted", "none", "good", "best"),
            nsmall = 300,
            renormalise = TRUE, normpower = 2,
            update = TRUE, leaveoneout = TRUE,
            reciplambda = NULL, sigma = NULL, adjust.sigma = 1, varcov = NULL,
            gref = NULL, tau = 0, fast = TRUE, var.approx = FALSE,
            domain = NULL, ratio = FALSE, close = NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{\dots}{
    Arguments passed to \code{\link[stats]{density.default}} or
    to \code{\link[spatstat.univar]{densityBC}} controlling the kernel smoothing.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
  }
  \item{adaptive}{
    Logical value specifying whether to use adaptive kernel smoothing
    (\code{adaptive=TRUE}) or fixed-bandwidth kernel smoothing
    (\code{adaptive=FALSE}, the default).
  }
  \item{kernel}{
    Choice of smoothing kernel, passed to \code{\link{density.default}}.
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
  \item{nsmall}{
    Optional. Integer. The maximum number of data points
    for which the default value of \code{zerocor} will be
    \code{"JonesFoster"}. 
  }
  \item{renormalise}{
    Logical. Whether to renormalise the estimate. See Details.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See Details.
  }
  \item{update}{
    Logical. If \code{lambda} is a fitted model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
    and \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link[spatstat.model]{update.ppm}} or \code{\link[spatstat.model]{update.kppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link[spatstat.model]{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
  \item{reciplambda}{
    Alternative to \code{lambda}.
    Values of the estimated \emph{reciprocal} \eqn{1/\lambda}{1/lambda}
    of the intensity function.
    Either a vector giving the reciprocal intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    reciprocal intensity values at all locations,
    or a \code{function(x,y)} which can be evaluated to give the
    reciprocal intensity value at any location.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
  }
  \item{adjust.sigma}{
    Numeric value. \code{sigma} will be multiplied by this value.
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
\details{
  The inhomogeneous pair correlation function \eqn{g_{\rm inhom}(r)}{ginhom(r)}
  is a summary of the dependence between points in a spatial point
  process that does not have a uniform density of points.

  The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda(x) lambda(y) g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda(x) * lambda(y) * g(r) dx dy
  }
  where \eqn{\lambda}{lambda} is the intensity function
  of the point process.
  For a Poisson point process with intensity function
  \eqn{\lambda}{lambda}, this probability is 
  \eqn{p(r) = \lambda(x) \lambda(y)}{p(r) = lambda(x) * lambda(y)}
  so \eqn{g_{\rm inhom}(r) = 1}{ginhom(r) = 1}.

  The inhomogeneous pair correlation function 
  is related to the inhomogeneous \eqn{K} function through
  \deqn{
    g_{\rm inhom}(r) = \frac{K'_{\rm inhom}(r)}{2\pi r}
  }{
    ginhom(r) = Kinhom'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'_{\rm inhom}(r)}{Kinhom'(r)}
  is the derivative of \eqn{K_{\rm inhom}(r)}{Kinhom(r)}, the
  inhomogeneous \eqn{K} function. See \code{\link{Kinhom}} for information
  about \eqn{K_{\rm inhom}(r)}{Kinhom(r)}.

  The command \code{pcfinhom} estimates the inhomogeneous
  pair correlation using a modified version of
  the algorithm in \code{\link{pcf}}. In this modified version,
  the contribution from each pair of points \eqn{X[i], X[j]} is
  weighted by
  \eqn{1/(\lambda(X[i]) \lambda(X[j]))}{1/(lambda(X[i]) lambda(X[j]))}.
  The arguments \code{divisor},
  \code{correction}, \code{zerocor} and \code{nsmall}
  are interpreted as described in the help file for \code{\link{pcf}}.
  
  If \code{renormalise=TRUE} (the default), then the estimates 
  are multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
  \eqn{
    c = \mbox{area}(W)/\sum (1/\lambda(x_i)).
  }{
    c = area(W)/sum[i] (1/lambda(x[i])).
  }
  This rescaling reduces the variability and bias of the estimate
  in small samples and in cases of very strong inhomogeneity.
  The default value of \code{normpower} is 2,
  which corresponds to rescaling
  the \code{lambda} values so that
  \eqn{
    \sum (1/\lambda(x_i)) = \mbox{area}(W).
  }{
    sum[i] (1/lambda(x[i])) = area(W).
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
\references{
  \smoothpcfpaper
}
\author{
  \spatstatAuthorsComma, \martinH and \tilman.
}
\seealso{
  \code{\link{pcf}},
  \code{\link{bw.bdh}},
  \code{\link{bw.pcfinhom}}
}
\examples{
  g <- pcfinhom(japanesepines, divisor="a")
}
\keyword{spatial}
\keyword{nonparametric}
