\name{Kmulti.inhom}
\alias{Kmulti.inhom}
\title{
  Inhomogeneous Marked K-Function
}
\description{
  For a marked point pattern, 
  estimate the inhomogeneous version of the multitype \eqn{K} function
  which counts the expected number of points of subset \eqn{J}
  within a given distance from a typical point in subset \code{I},
  adjusted for spatially varying intensity.
}
\usage{
  Kmulti.inhom(X, I, J, lambdaI=NULL, lambdaJ=NULL,
          \dots,
          r=NULL, breaks=NULL, rmax=NULL,
          correction=c("border", "isotropic", "Ripley", "translate"),
          lambdaIJ=NULL,
          sigma=NULL, varcov=NULL,
          lambdaX=NULL, update=TRUE, leaveoneout=TRUE,
          Ilevels=NULL, Jlevels=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous multitype \eqn{K} function
    \eqn{K_{IJ}(r)}{KIJ(r)} will be computed.
    It must be a marked point pattern.
    See under Details.
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured. See Details.
  }
  \item{J}{Subset index specifying the points in \code{X} to which
    distances are measured. See Details.
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
  \item{\dots}{Ignored.}
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the multitype \eqn{K} function
    \eqn{K_{IJ}(r)}{KIJ(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
    A single numeric value. There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{lambdaIJ}{
    Optional. A matrix containing estimates of
    the product of the intensities \code{lambdaI} and \code{lambdaJ}
    for each pair of points, the first point belonging to subset
    \code{I} and the second point to subset \code{J}.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
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
  \item{Ilevels}{
    Character vector containing the types of points
    (levels of \code{marks(X)}) which
    comprise the subset \code{I}.
    Required when \code{lambdaX} is given
    or when \code{lambdaI} is a multitype point process model.
  }
  \item{Jlevels}{
    Character vector containing the types of points
    which comprise the subset \code{J}.
    Required when \code{lambdaX} is given
    or when \code{lambdaJ} is a multitype point process model.
  }
}
\value{
An object of class \code{"fv"} (see \code{\link{fv.object}}).

Essentially a data frame containing numeric columns 
\item{r}{the values of the argument \eqn{r} 
at which the function \eqn{K_{IJ}(r)}{KIJ(r)} has been  estimated
}
\item{theo}{the theoretical value of  \eqn{K_{IJ}(r)}{KIJ(r)}
for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
}
together with a column or columns named 
\code{"border"}, \code{"bord.modif"},
\code{"iso"} and/or \code{"trans"},
according to the selected edge corrections. These columns contain
estimates of the function \eqn{K_{IJ}(r)}{KIJ(r)}
obtained by the edge corrections named.
}
\details{
  The function \code{Kmulti.inhom}
  is the counterpart, for spatially-inhomogeneous marked point patterns,
  of the multitype \eqn{K} function \code{\link{Kmulti}}.

  Suppose \eqn{X} is a marked point process, with marks of any kind.
  Suppose \eqn{X_I}{X[I]}, \eqn{X_J}{X[J]} are two sub-processes, possibly
  overlapping. Typically \eqn{X_I}{X[I]} would consist of those points
  of \eqn{X} whose marks lie in a specified range of mark values,
  and similarly for \eqn{X_J}{X[J]}. Suppose that
  \eqn{\lambda_I(u)}{lambdaI(u)}, \eqn{\lambda_J(u)}{lambdaJ(u)} are the
  spatially-varying intensity functions of \eqn{X_I}{X[I]} and
  \eqn{X_J}{X[J]} respectively. Consider all the pairs of points
  \eqn{(u,v)} in the point process \eqn{X} such that the first point
  \eqn{u} belongs to \eqn{X_I}{X[I]}, the second point \eqn{v}
  belongs to \eqn{X_J}{X[J]}, and the distance between \eqn{u} and \eqn{v}
  is less than a specified distance \eqn{r}. Give this pair \eqn{(u,v)}
  the numerical weight
  \eqn{1/(\lambda_I(u)\lambda_J(u))}{1/(lambdaI(u) lambdaJ(u))}.
  Calculate the sum of these weights over all pairs of points as
  described. This sum (after appropriate edge-correction and
  normalisation) is the estimated inhomogeneous multitype \eqn{K} function.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link[spatstat.geom]{as.ppp}}.

  The arguments \code{I} and \code{J} specify two subsets of the
  point pattern. They may be any type of subset indices, for example,
  logical vectors of length equal to \code{npoints(X)},
  or integer vectors with entries in the range 1 to
  \code{npoints(X)}, or negative integer vectors.

  Alternatively, \code{I} and \code{J} may be \bold{functions}
  that will be applied to the point pattern \code{X} to obtain
  index vectors. If \code{I} is a function, then evaluating
  \code{I(X)} should yield a valid subset index. This option
  is useful when generating simulation envelopes using
  \code{\link{envelope}}.

  The argument \code{lambdaI} supplies the values
  of the intensity of the sub-process identified by index \code{I}.
  It may be either
  \describe{
    \item{a pixel image}{(object of class \code{"im"}) which
      gives the values of the intensity of \code{X[I]}
      at all locations in the window containing \code{X};
    }
    \item{a numeric vector}{containing the values of the
      intensity of \code{X[I]} evaluated only
      at the data points of \code{X[I]}. The length of this vector
      must equal the number of points in \code{X[I]}.
    }
    \item{a function}{
      of the form \code{function(x,y)}
      which can be evaluated to give values of the intensity at
      any locations.
    }
    \item{A fitted multitype point process model}{
      (class \code{"ppm"})
      that was fitted to a point pattern. If
      \code{update=TRUE}, this model will be re-fitted to the data
      pattern \code{X}. Then the intensity of this model for
      the subset \code{I} will be computed (using the information
      \code{Ilevels}).
    }
    \item{A fitted unmarked point process model}{
      (class \code{"ppm"}, \code{"kppm"}, \code{"dppm"} or \code{"slrm"})
      that was fitted to an \emph{unmarked} point pattern. If
      \code{update=TRUE}, this model will be re-fitted to the subset of
      the data, \code{X[I]}. Then the intensity of this fitted model
      will be computed.
    }
    \item{omitted:}{
      if \code{lambdaI} is omitted then it will be estimated
      using a leave-one-out kernel smoother. 
    }
  }
  If \code{lambdaI} is omitted, then it will be estimated using
  a `leave-one-out' kernel smoother, as described in Baddeley,
  \Moller and Waagepetersen (2000).  The estimate of \code{lambdaI} for a given
  point is computed by removing the point from the
  point pattern, applying kernel smoothing to the remaining points using
  \code{\link{density.ppp}}, and evaluating the smoothed intensity
  at the point in question. The smoothing kernel bandwidth is controlled
  by the arguments \code{sigma} and \code{varcov}, which are passed to
  \code{\link{density.ppp}} along with any extra arguments.

  Similarly \code{lambdaJ} supplies the values
  of the intensity of the sub-process identified by index \code{J}.

  Alternatively if the argument \code{lambdaX} is given,
  then it should
  provide data giving the intensities of both subsets \code{I} and
  \code{J}. It may be
  \itemize{
    \item A list of pixel images, one for each possible type of point
    \item A \code{function(x,y,mark)} which can be evaluated to
    give the intensity of each type of point at any location
    \item A matrix with one row for each data point in \code{X}
    and one column for each possible type of point, giving the estimated
    intensity of each type of point at each data location
    \item A fitted point process model (class \code{"ppm"})
    that was fitted to a multitype point pattern. If
    \code{update=TRUE}, this model will be re-fitted to the data
    pattern \code{X}. 
  }
  The arguments \code{lambdaI}, \code{lambdaJ} will be ignored
  if \code{lambdaX} is given.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{IJ}(r)}{KIJ(r)} should be evaluated. 
  It is also used to determine the breakpoints
  (in the sense of \code{\link[graphics]{hist}})
  for the computation of histograms of distances.

  First-time users would be strongly advised not to specify \code{r}.
  However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, 
  and \code{max(r)} must be larger than the radius of the largest disc 
  contained in the window. 

  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kinhom}}.
  The edge corrections implemented here are
  \describe{
    \item{border}{the border method or
      ``reduced sample'' estimator (see Ripley, 1988). This is
      the least efficient (statistically) and the fastest to compute.
      It can be computed for a window of arbitrary shape.
    }
    \item{isotropic/Ripley}{Ripley's isotropic correction
      (see Ripley, 1988; Ohser, 1983).
      This is currently implemented only for rectangular windows.
    }
    \item{translate}{Translation correction (Ohser, 1983).
      Implemented for all window geometries.
    }
  }
  The pair correlation function \code{\link{pcf}} can also be applied to the
  result of \code{Kmulti.inhom}.
}
\references{
  Baddeley, A., \Moller, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.
}
\seealso{
 \code{\link{Kmulti}},
 \code{\link{Kdot.inhom}},
 \code{\link{Kcross.inhom}},
 \code{\link{pcf}}
}
\examples{
    # Finnish Pines data: marked by diameter and height
    plot(finpines, which.marks="height")
    II <- (marks(finpines)$height <= 2)
    JJ <- (marks(finpines)$height > 3)
    K <- Kmulti.inhom(finpines, II, JJ)
    plot(K)
    # functions determining subsets
    f1 <- function(X) { marks(X)$height <= 2 }
    f2 <- function(X) { marks(X)$height > 3 }
    K <- Kmulti.inhom(finpines, f1, f2)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}



