\name{pcfmulti.inhom}
\alias{pcfmulti.inhom}
\title{
  Inhomogeneous Multitype Pair Correlation Function
}
\description{
  Estimates the inhomogeneous multitype pair correlation function
  for a multitype point pattern.
}
\usage{
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, \dots,
               lambdaX=NULL,
               r = NULL, breaks = NULL, rmax = NULL,
               adaptive = FALSE,
               kernel = "epanechnikov", bw = NULL, h = NULL, bw.args = list(),
               stoyan = 0.15, adjust.bw = 1,
               correction = c("translate", "Ripley"),
               divisor = c("a", "r", "d", "t"),
               zerocor = c("convolution", "reflection",
                         "bdrykern", "JonesFoster", "weighted",
                         "none", "good", "best"),
               nsmall = 300,
               gref = NULL, tau = 0,
               sigma = NULL, adjust.sigma = 1, varcov = NULL,
               update = TRUE, leaveoneout = TRUE,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J",
               IJexclusive = FALSE, Ilevels=NULL, Jlevels=NULL,
               close = NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous
    cross-type pair correlation function
    \eqn{g_{ij}(r)}{g[i,j](r)}
    will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). 
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured. 
  }
  \item{J}{Subset index specifying the points in \code{X}
    to which distances are measured. 
  }
 \item{lambdaI}{
    Optional.
    Values of the estimated intensity function of the points
    belonging to subset \code{I}.
    (Ignored if \code{lambdaX} is given.)
    Either a numeric vector giving the intensity values
    at the data points of subset \code{I},
    a pixel image (object of class \code{"im"}) giving the
    intensity values of points of subset \code{I} at all locations,
    a \code{function(x,y)} which
    can be evaluated to give the intensity value of points of subset \code{I}
    at any location, or a fitted point process model
    (class \code{"ppm"}, \code{"kppm"}, \code{"dppm"} or \code{"slrm"})
    which could be used to predict the
    intensity values of points of subset \code{J} at all locations.
    If \code{lambdaI} is an unmarked point process model, it is assumed
    to have been fitted to \code{X[I]}; if it is a multitype point
    process model, it is assumed to have been fitted to \code{X},
    and the argument \code{Ilevels} is required.
  }
  \item{lambdaJ}{
    Optional.
    Values of the estimated intensity function of the points
    belonging to subset \code{J}.
    A numeric vector, pixel image, \code{function(x,y)},
    or fitted point process model.
    Ignored if \code{lambdaX} is given.
  }
  \item{\dots}{
    Ignored.
  }
  \item{lambdaX}{
    Optional.
    Alternative to \code{lambdaI} and \code{lambdaJ}.
    Data which can be used to calculate the intensity functions
    of both the subsets \code{I} and \code{J}.
    A list of pixel images (one image for each possible type of point),
    a \code{function(x,y,mark)} giving the intensity at each location
    for each possible type of point, a numeric matrix with one row for
    each point in \code{X} and one column for each possible type of
    point, or a fitted multitype point process model (class
    \code{"ppm"}). The arguments \code{Ilevels} and \code{Jlevels} are
    required in this case.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{breaks}{Internal use only.}
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
  \item{adjust.bw}{
    Numerical adjustment factor for the bandwidth.
    The bandwidth actually used is \code{bw * adjust.bw}.
    This makes it easy to specify choices like
    \sQuote{half the selected bandwidth}.
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
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
  }
  \item{adjust.sigma}{
    Numeric value. \code{sigma} will be multiplied by this value.
  }
  \item{update}{
    Logical value indicating what to do when
    \code{lambdaI}, \code{lambdaJ} or \code{lambdaX}
    is a fitted point process model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
    If \code{update=TRUE} (the default),
    the model will first be updated (refitted to the data)
    before the fitted intensity \code{lambdaI} or \code{lambdaJ}
    at the data points is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting the model.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link[spatstat.model]{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
  \item{Iname,Jname}{
    Optional. Character strings describing the members of
    the subsets \code{I} and \code{J}.
  }
  \item{IJexclusive}{
    Logical value indicating whether the subsets \code{I} and \code{J}
    are guaranteed to be mutually exclusive.
  }
  \item{Ilevels}{
    Character vector containing the types of points
    (levels of \code{marks(X)}) which
    comprise the subset \code{I}.
    Required when \code{lambdaX} is given
    or when \code{lambdaI} is a multitype point process model.
  }
  \item{Jlevels}{
    Character vector containing the types of points
    (levels of \code{marks(X)}) which
    comprise the subset \code{J}.
    Required when \code{lambdaX} is given
    or when \code{lambdaJ} is a multitype point process model.
  }
  \item{close}{
    Advanced use only.
    Precomputed data obtained from \code{\link[spatstat.geom]{crosspairs}}.
  }
}
\details{
  This is a generalisation of \code{\link{pcfcross.inhom}}
  to arbitrary collections of points.

  The algorithm measures the distance from each data point
  in subset \code{I} to each data point in subset \code{J},
  excluding identical pairs of points. The distances are
  weighted by the reciprocal intensity (as described below), then
  kernel-smoothed and renormalised to form a pair correlation
  function.

  The contribution from each pair of data points \eqn{u,v} is weighted by
  \eqn{1/(\lambda_I(u) \lambda_J(v))}{1/(lambdaI(u) * lambdaJ(v))}
  where \eqn{\lambda_I(u)}{lambdaI(u)} is the estimated intensity
  of points belonging to subset \code{I} at location \eqn{u},
  and similarly \eqn{\lambda_J(u)}{lambdaJ(u)} is the estimated intensity
  of points belonging to subset \code{J} at location \eqn{v}.

  The estimated intensities \eqn{\lambda_I(u)}{lambdaI(u)}
  and \eqn{\lambda_J(u)}{lambdaJ(u)} are determined as follows.
  \itemize{
    \item
    If \code{lambdaX} is given, then it should
    provide data giving the intensities of both subsets \code{I} and
    \code{J}. It may be
    \itemize{
      \item A list of pixel images, one for each possible type of point
      \item A \code{function(x,y,mark)} which can be evaluated to
      give the intensity of each type of point at any location
      \item A matrix with one row for each data point in \code{X}
      and one column for each possible type of point, giving the estimated
      intensity of each type of point at each data location
      \item A fitted point process model (class \code{"lppm"} or \code{"ppm"})
      that was fitted to a multitype point pattern. If
      \code{update=TRUE}, this model will be re-fitted to the data
      pattern \code{X}. Then the fitted intensities for subsets \code{I} and
      \code{J} will be computed (using the information \code{Ilevels}
      and \code{Jlevels}).
    }
    \item
    Otherwise if \code{lambdaI} is given, it may be
    \itemize{
      \item A pixel image
      \item A \code{function(x,y)}
      \item A numeric vector with one entry for each data point in \code{X[I]}
      \item A fitted point process model (class \code{"lppm"} or \code{"ppm"})
      that was fitted to a \emph{multitype} point pattern. If
      \code{update=TRUE}, this model will be re-fitted to the data
      pattern \code{X}. Then the intensity of this model for
      the subset \code{I} will be computed (using the information
      \code{Ilevels}).
      \item A fitted point process model
      (class \code{"lppm"}, \code{"ppm"}, \code{"kppm"}, \code{"dppm"} or \code{"slrm"})
      that was fitted to an \emph{unmarked} point pattern. If
      \code{update=TRUE}, this model will be re-fitted to the relevant subset of
      the data, \code{unmark(X[I])}. Then the intensity of this fitted model
      will be computed.
    }
    \item
    Otherwise if \code{lambda} is missing or \code{NULL},
    then a kernel estimate of the intensity
    \eqn{\lambda_I(u)}{lambdaI(u)} 
    will be computed by applying \code{\link{density.ppp}}
    to \code{X[I]} using the arguments
    \code{sigma}, \code{adjust.sigma}, \code{varcov}
    and \code{leaveoneout}.
    \item
    Similarly for \code{lambdaJ}.
  }
  
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
  \code{\link{pcfcross.inhom}}, 
  \code{\link{pcf.ppp}}.
}
\author{
  \spatstatAuthorsComma, \tilman and \martinH.
}
\examples{
  adult <- (marks(longleaf) >= 30)
  juvenile <- !adult
  p <- pcfmulti.inhom(longleaf, adult, juvenile)
}
\keyword{spatial}
\keyword{nonparametric}
