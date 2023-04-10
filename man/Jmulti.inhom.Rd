\name{Jmulti.inhom}
\alias{Jmulti.inhom}
\title{
  Inhomogeneous Marked J-Function
}
\description{
  For a marked point pattern, 
  estimate the inhomogeneous version of the multitype \eqn{J} function.
}
\usage{
Jmulti.inhom(X, I, J,
             lambda = NULL, lambdaI = NULL, lambdaJ = NULL,
             lambdamin = NULL,
             \dots,
             r = NULL,
             ReferenceMeasureMarkSetI = NULL,
             ratio = FALSE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous multitype \eqn{J} function
    \eqn{J_{IJ}(r)}{J[IJ](r)} will be computed.
    It must be a marked point pattern.
    See under Details.
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured, for the inhomogeneous \eqn{G} function.
  }
  \item{J}{Subset index specifying the points in \code{X} to which
    distances are measured, for the inhomogeneous \eqn{G} and \eqn{F} functions.
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}) or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
 \item{lambdaI}{
    Optional.
    Values of the estimated intensity of the sub-process \code{X[I]}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the points in \code{X[I]},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location,
  }
 \item{lambdaJ}{
    Optional.
    Values of the estimated intensity of the sub-process \code{X[J]}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the points in \code{X[J]},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lambdamin}{
    Optional. The minimum possible value of the intensity
    over the spatial domain. A positive numerical value.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link[spatstat.geom]{as.mask}} to control
    the pixel resolution for the computation.
  }
  \item{r}{
    vector of values for the argument \eqn{r} at which
    the inhomogeneous \eqn{K} function
    should be evaluated. Not normally given by the user;
    there is a sensible default.
  }
  \item{ReferenceMeasureMarkSetI}{
    Optional. The total measure of the mark set. A positive number.
  }
  \item{ratio}{
    Logical value indicating whether to save ratio information.
  }
}
\details{
  This function is the counterpart of \code{\link{Jmulti}}
  for inhomogeneous patterns. It is computed by evaluating the
  inhomogeneous \eqn{G} function \code{\link{GmultiInhom}} 
  and the inhomogeneous \eqn{F} function \code{\link{FmultiInhom}}
  and computing the ratio \eqn{J = (1-G)/(1-F)}.
}
\value{
  Object of class \code{"fv"} containing the estimate of the
  inhomogeneous multitype \eqn{J} function.
}
\references{
  Cronie, O. and Van Lieshout, M.N.M. (2015)
  Summary statistics for inhomogeneous marked point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  DOI: 10.1007/s10463-015-0515-z
}
\author{
  Jonatan \Gonzalez and \adrian.
}
\seealso{
  \code{\link{Jcross.inhom}}, \code{\link{Jdot.inhom}}
  for special cases.
  
  \code{\link{GmultiInhom}}, \code{\link{FmultiInhom}},
  \code{\link{Jmulti}}.
}
\examples{
  X <- rescale(amacrine)
  I <- (marks(X) == "on")
  J <- (marks(X) == "off")
  if(interactive() && require(spatstat.model)) {
    ## how to do it normally
    mod <- ppm(X ~ marks * x)
    lam <- fitted(mod, dataonly=TRUE)
    lmin <- min(predict(mod)[["off"]]) * 0.9
    dd <- NULL
  } else {
    ## for package testing
    lam <- intensity(X)[as.integer(marks(X))]
    lmin <- intensity(X)[2] * 0.9
    dd <- 32
  }
  JM <- Jmulti.inhom(X, I, J, lambda=lam, lambdamin=lmin, dimyx=dd)
}
\keyword{spatial}
\keyword{nonparametric}
