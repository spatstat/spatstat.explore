\name{pcfcross.inhom}
\alias{pcfcross.inhom}
\title{
  Inhomogeneous Multitype Pair Correlation Function (Cross-Type)
}
\description{
  Estimates the inhomogeneous cross-type pair correlation function
  for a multitype point pattern.
}
\usage{
pcfcross.inhom(X, i, j, lambdaI = NULL, lambdaJ = NULL, \dots) 
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
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be
    converted to a character string).
    Defaults to the second level of \code{marks(X)}.
  }
  \item{lambdaI}{
    Optional.
    Values of the estimated intensity function of the points of type \code{i}.
    Either a numeric vector giving the intensity values
    at the data points of type \code{i},
    a pixel image (object of class \code{"im"}) giving the
    intensity values of points of type \code{i} at all locations,
    a \code{function(x,y)} which
    can be evaluated to give the intensity value of points of type
    \code{i} at any location, a fitted multitype point process model
    (class \code{"ppm"}) 
    which could be used to predict the
    intensity values of points of \emph{each} type at any location,
    or a fitted unmarked point process model
    (class \code{"ppm"}, \code{"kppm"}, \code{"dppm"} or \code{"slrm"})
    which could be used to predict the
    intensity values of points of type \code{i} only at any location.
  }
  \item{lambdaJ}{
    Optional.
    Values of the estimated intensity function of the points of type \code{j}.
    Numeric vector, pixel image, \code{function(x,y)} or fitted point
    process model.
  }
  \item{\dots}{
    Arguments passed to \code{\link{pcfmulti.inhom}} to control the
    computation. 
  }
}
\details{
  The inhomogeneous cross-type pair correlation function
  \eqn{g_{ij}(r)}{g[i,j](r)}
  is a summary of the dependence between two types of points in a
  multitype spatial point process that does not have a uniform
  density of points.

  The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points, of types \eqn{i} and \eqn{j}
  respectively, at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda_i(x) lambda_j(y) g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda[i](x) * lambda[j](y) * g(r) dx dy
  }
  where \eqn{\lambda_i}{lambda[i]} is the intensity function
  of the process of points of type \eqn{i}.
  For a multitype Poisson point process,
  this probability is 
  \eqn{p(r) = \lambda_i(x) \lambda_j(y)}{p(r) = lambda[i](x) * lambda[j](y)}
  so   \eqn{g_{ij}(r) = 1}{g[i,j](r) = 1}.

  The command \code{pcfcross.inhom} estimates the inhomogeneous
  pair correlation using a modified version of
  the algorithm in \code{\link{pcf.ppp}}.
  The arguments \code{bw} and \code{adjust.bw} control the
  degree of one-dimensional smoothing of the estimate of pair correlation.

  If the arguments \code{lambdaI} and/or \code{lambdaJ} are missing or
  null, they will be estimated from \code{X} by spatial kernel smoothing
  using a leave-one-out estimator, computed by \code{\link{density.ppp}}.
  The arguments \code{sigma}, \code{varcov}
  and \code{adjust.sigma} control the degree of spatial smoothing.
}
\value{
  A function value table (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{
    the vector of values of the argument \eqn{r} 
    at which the inhomogeneous cross-type pair correlation function
    \eqn{g_{ij}(r)}{g[i,j](r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g_{ij}(r)}{g[i,j](r)}
    for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g_{ij}(r)}{g[i,j](r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g_{ij}(r)}{g[i,j](r)}
    estimated by Ripley isotropic correction
  }
  as required.
}
\seealso{
  \code{\link{pcf.ppp}}, 
  \code{\link{pcfinhom}}, 
  \code{\link{pcfcross}},
  \code{\link{pcfdot.inhom}}
}
\examples{
  plot(pcfcross.inhom(amacrine, "on", "off", stoyan=0.1),
       legendpos="bottom")
}
\author{\spatstatAuthorsComma, \tilman and \martinH.}
\keyword{spatial}
\keyword{nonparametric}
