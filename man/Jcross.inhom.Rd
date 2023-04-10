\name{Jcross.inhom}
\alias{Jcross.inhom}
\title{
  Inhomogeneous Multitype J function (i-to-j)
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous multitype \eqn{J} function 
  summarising the interpoint dependence between
  points of type \eqn{i} and of type \eqn{j}.
}
\usage{
  Jcross.inhom(X, i, j,
               lambda = NULL, lambdaI = NULL, lambdaJ = NULL,
               lambdamin = NULL,
               \dots,
               r = NULL, ReferenceMeasureMarkSetI = NULL, ratio = FALSE)
}
\arguments{
  \item{X}{
    The observed point pattern, 
    from which an estimate of the multitype \eqn{J} function
    \eqn{J_{ij}(r)}{Jij(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the second level of \code{marks(X)}.
  }
 \item{lambda}{
    Optional.
    Values of the estimated intensity of the point process.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the points in \code{X},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lambdaI}{
    Optional.
    Values of the estimated intensity of the sub-process of
    points of type \code{i}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the type \code{i} points in \code{X},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lambdaJ}{
    Optional.
    Values of the the estimated intensity of the sub-process of
    points of type \code{j}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the type \code{j} points in \code{X},
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
    the inhomogeneous \eqn{J} function
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
  This function is the counterpart of \code{\link{Jcross}}
  for inhomogeneous patterns. It is computed as a special case
  of \code{\link{Jmulti.inhom}}.
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
  \code{\link{Jdot.inhom}}, \code{\link{Jmulti.inhom}},
  \code{\link{Jcross}}.
}
\examples{
  X <- rescale(amacrine)
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
  JC <- Jcross.inhom(X, "on", "off", lambda=lam, lambdamin=lmin, dimyx=dd)
}
