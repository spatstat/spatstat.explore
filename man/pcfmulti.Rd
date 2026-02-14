\name{pcfmulti}
\alias{pcfmulti}
\title{
  Marked pair correlation function
}
\description{
For a marked point pattern, 
estimate the multitype pair correlation function
using kernel methods.
}
\usage{
   pcfmulti(X, I, J, ..., r = NULL, rmax=NULL,
            adaptive=FALSE,
            kernel = "epanechnikov", bw = NULL, h=NULL,
            bw.args=list(), stoyan = 0.15, adjust=1,
            correction = c("translate", "Ripley"),
            divisor=c("a", "r", "d", "t"),
            zerocor=c("convolution", "reflection", "bdrykern",
                      "JonesFoster", "weighted", "none",
                      "good", "best"),
            nsmall = 300,
            gref=NULL,
            tau = 0,
            Iname = "points satisfying condition I",
            Jname = "points satisfying condition J",
            IJexclusive=FALSE,
            ratio = FALSE,
            close=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype pair correlation function
    \eqn{g_{IJ}(r)}{g[I,J](r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). 
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured. 
  }
  \item{J}{Subset index specifying the points in \code{X}
    to which distances are measured. 
  }
  \item{\dots}{
    Ignored.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
    A single numeric value. There is a sensible default.
  }
  \item{adaptive}{
    Logical value specifying whether to use adaptive kernel smoothing
    (\code{adaptive=TRUE}) or fixed-bandwidth kernel smoothing
    (\code{adaptive=FALSE}, the default).
  }
  \item{kernel}{
    Choice of smoothing kernel,
    passed to \code{\link{density.default}}.
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
    Coefficient for default bandwidth rule.
  }
  \item{adjust}{
     Numerical adjustment factor for the bandwidth.
     The bandwidth actually used is \code{adjust * bw}.
     This makes it easy to specify choices like \sQuote{half the
     selected bandwidth}.
  }
  \item{correction}{
    String (partially matched) specifying the choice or choices
    of spatial edge correction. Options include \code{"translate"} for
    the translation correction, \code{"isotropic"} or \code{"Ripley"}
    for Ripley's isotropic correction, and \code{"none"} for no edge correction.
  }
  \item{divisor}{
    String specifying the choice of estimator.
    See \code{\link{pcf.ppp}}.
  }
  \item{zerocor}{
    String (partially matched) specifying a correction for the boundary effect
    bias at \eqn{r=0}. Possible values are
    \code{"none"}, \code{"weighted"}, \code{"convolution"},
    \code{"reflection"}, \code{"bdrykern"} and \code{"JonesFoster"}.
    See \code{\link{pcf.ppp}}.
  }
  \item{nsmall}{
    Optional. Integer. The maximum number of data points
    for which the default value of \code{zerocor} will be
    \code{"JonesFoster"}. 
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
  \item{Iname,Jname}{
    Optional. Character strings describing the members of
    the subsets \code{I} and \code{J}.
  }
  \item{IJexclusive}{
    Logical value indicating whether the subsets \code{I} and \code{J}
    are guaranteed to be mutually exclusive.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
  \item{close}{
    Advanced use only. Precomputed data obtained from \code{\link[spatstat.geom]{crosspairs}}.
  }
}
\details{
  This is a generalisation of \code{\link{pcfcross}}
  to arbitrary collections of points.

  The algorithm measures the distance from each data point
  in subset \code{I} to each data point in subset \code{J},
  excluding identical pairs of points. The distances are
  kernel-smoothed and renormalised to form a pair correlation
  function.

  The smoothing algorithm is a multitype version of the
  smoothing algorithm in \code{\link{pcf.ppp}}.

  See \code{\link{pcf.ppp}} for detailed documentation of the arguments
  \code{correction}, \code{divisor}, \code{zerocor},
  and other smoothing arguments.
}
\value{
  An object of class \code{"fv"}.
}
\seealso{
  \code{\link{pcfcross}}, 
  \code{\link{pcfdot}}, 
  \code{\link{pcf.ppp}}.
}
\examples{
  adult <- (marks(longleaf) >= 30)
  juvenile <- !adult
  p <- pcfmulti(longleaf, adult, juvenile)
}
\author{
  \spatstatAuthorsComma, \tilman and \martinH.
}
\keyword{spatial}
\keyword{nonparametric}

