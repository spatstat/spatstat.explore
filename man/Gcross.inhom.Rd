\name{Gcross.inhom}
\alias{Gcross.inhom}
\title{
  Inhomogeneous Multitype G Cross Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the cross \eqn{G} function,
  which is the distribution of the distance
  from a point of type \eqn{i}
  to the nearest point of type \eqn{j},
  adjusted for spatially varying intensity.
}
\usage{
   Gcross.inhom(X, i, j,
              lambda = NULL, lambdaI = NULL, lambdaJ = NULL,
              lambdamin = NULL,
              \dots,
              r = NULL,
              ReferenceMeasureMarkSetI = NULL,
              ratio = FALSE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{G} function
    \eqn{G_{ij}(r)}{Gij(r)} will be computed.
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
    the inhomogeneous \eqn{G} function
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
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}})
  containing estimates of the inhomogeneous cross type \eqn{G} function.
}
\details{
  This is a generalisation of the function \code{\link{Gcross}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Ginhom}}.
  
  The argument \code{lambdaI} supplies the values
  of the intensity of the sub-process of points of type \code{i}.
  It may be either
  \describe{
    \item{a pixel image}{(object of class \code{"im"}) which
      gives the values of the type \code{i} intensity
      at all locations in the window containing \code{X};
    }
    \item{a numeric vector}{containing the values of the
      type \code{i} intensity evaluated only
      at the data points of type \code{i}. The length of this vector
      must equal the number of type \code{i} points in \code{X}.
    }
    \item{a function}{
      of the form \code{function(x,y)}
      which can be evaluated to give values of the intensity at
      any locations.
    }
     \item{a fitted point process model}{
      (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
      whose fitted \emph{trend} can be used as the fitted intensity.
      (If \code{update=TRUE} the model will first be refitted to the
      data \code{X} before the trend is computed.)
    }
   \item{omitted:}{
      if \code{lambdaI} is omitted then it will be estimated
      using a leave-one-out kernel smoother. 
    }
  }
  If \code{lambdaI} is omitted, then it will be estimated using
  a `leave-one-out' kernel smoother.

  Similarly the argument \code{lambdaJ} should contain
  estimated values of the intensity of the points of type \eqn{j}.
  It may be either a pixel image, a numeric vector of length equal
  to the number of points in \code{X}, a function, or omitted.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{G_{ij}(r)}{Gij(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must not exceed the radius of the
  largest disc contained in the window.
}
\references{
  Cronie, O. and Van Lieshout, M.N.M. (2015)
  Summary statistics for inhomogeneous marked point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  DOI: 10.1007/s10463-015-0515-z
}
\section{Warnings}{
  The argument \code{i} is interpreted as
  a level of the factor \code{X$marks}. It is converted to a character
  string if it is not already a character string.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.
}
\seealso{
 \code{\link{Gcross}},
 \code{\link{Ginhom}},
 \code{\link{Gcross.inhom}},
 \code{\link{Gmulti.inhom}}.
}
\examples{
  X <- rescale(amacrine)
  if(interactive() && require(spatstat.model)) {
    ## how to do it normally
    mod <- ppm(X ~ marks * x)
    lam <- fitted(mod, dataonly=TRUE)
    lmin <- min(predict(mod)[["off"]]) * 0.9
  } else {
    ## for package testing 
    lam <- intensity(X)[as.integer(marks(X))]
    lmin <- intensity(X)[2] * 0.9
  }
  GC <- Gcross.inhom(X, "on", "off", lambda=lam, lambdamin=lmin)
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{nonparametric}

