\name{relriskHeat}
\alias{relriskHeat}
\alias{relriskHeat.ppp}
\title{
  Diffusion Estimate of Conditional Probabilities
}
\description{
  Computes the conditional probability estimator of relative risk
  based on a multitype point pattern using the diffusion estimate
  of the type-specific intensities.
}
\usage{
relriskHeat(X, \dots)

\method{relriskHeat}{ppp}(X, \dots, sigmaX=NULL, weights=NULL)
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{densityHeat}}
    controlling the estimation of each marginal intensity.
  }
  \item{sigmaX}{
    Optional.
    Numeric vector of bandwidths, one associated with each data point in
    \code{X}. 
  }
  \item{weights}{
    Optional numeric vector of weights associated with each point of
    \code{X}.
  }
}
\details{
  The function \code{relriskHeat} is generic. This file documents the
  method \code{relriskHeat.ppp} for spatial point patterns (objects of
  class \code{"ppp"}).
  
  This function estimates the spatially-varying conditional probability
  that a random point (given that it is present) will belong to
  a given type.

  The algorithm separates \code{X} into 
  the sub-patterns consisting of points of each type.
  It then applies \code{\link[spatstat.explore]{densityHeat}} to each sub-pattern,
  using the same bandwidth and smoothing regimen for each sub-pattern,
  as specified by the arguments \code{\dots}.

  If \code{weights} is specified, it should be a numeric vector
  of length equal to the number of points in \code{X}, so that
  \code{weights[i]} is the weight for data point \code{X[i]}.
  
  Similarly when performing lagged-arrival smoothing,
  the argument \code{sigmaX} must be a numeric vector of the same length
  as the number of points in \code{X}, and thus contain the
  point-specific bandwidths in the order corresponding to each of these
  points regardless of mark.
}
\value{
  A named list (of class \code{\link[spatstat.geom]{solist}})
  containing pixel \code{\link[spatstat.geom]{im}}ages,
  giving the estimated conditional probability surfaces for each type.
}
\seealso{
  \code{\link[spatstat.explore]{relrisk.ppp}} for the
  traditional convolution-based kernel estimator of
  conditional probability surfaces,
  and the function \code{risk} in the \pkg{sparr} package for the
  density-ratio-based estimator.
}
\references{
  Agarwal, N. and Aluru, N.R. (2010) 
  A data-driven stochastic collocation approach for
  uncertainty quantification in MEMS.
  \emph{International Journal for Numerical Methods in Engineering}
  \bold{83}, 575--597.

  Baddeley, A., Davies, T., Rakshit, S., Nair, G. and McSwiggan, G. (2022)
  Diffusion smoothing for spatial point patterns.
  \emph{Statistical Science} \bold{37}, 123--142.
  
  Barry, R.P. and McIntyre, J. (2011)
  Estimating animal densities and home range in regions with irregular
  boundaries and holes: a lattice-based alternative to the kernel
  density estimator. \emph{Ecological Modelling} \bold{222}, 1666--1672.

  Botev, Z.I. and Grotowski, J.F. and Kroese, D.P. (2010)
  Kernel density estimation via diffusion.
  \emph{Annals of Statistics} \bold{38}, 2916--2957.
}
\author{
  \adrian and \tilman.
}
\examples{
  ## bovine tuberculosis data
  X <- subset(btb, select=spoligotype)
  plot(X) 
  P <- relriskHeat(X,sigma=9)
  plot(P)
}
\keyword{spatial}
\keyword{smooth}
