\name{relriskAdaptiveKernel.ppp}
\alias{relriskAdaptiveKernel}
\alias{relriskAdaptiveKernel.ppp}
\title{
  Adaptive Estimate of Spatially-Varying Probability or Relative Risk
}
\description{
  Given a multitype point pattern, estimate the
  spatially-varying probability of each type of point, or the ratios of
  such probabilities, using kernel smoothing
  with adaptive bandwidths.
}
\usage{
relriskAdaptiveKernel(X, \dots)

\method{relriskAdaptiveKernel}{ppp}(X, bw = NULL, \dots,
      adjust = 1, at = c("pixels", "points"),
      weights = NULL, relative = FALSE, normalise = FALSE, 
      casecontrol = TRUE, control = 1, case)
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}
    which has factor valued marks) or a point pattern with several
    columns of factor-valued marks.
  }
  \item{bw}{
    Numeric vector of smoothing bandwidths for each point in \code{X},
    or a pixel image giving the smoothing bandwidth at each spatial
    location, or a spatial function of class \code{"funxy"} giving the
    smoothing bandwidth at each location.
    The default is to compute bandwidths using
    \code{\link[spatstat.explore]{bw.abram.ppp}}.
  }
  \item{adjust}{
    Optional. Adjustment factor for the bandwidths \code{bw}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{bw.abram.ppp}} to select the
    bandwidths, or passed to
    \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}} 
    to control the pixel resolution.
  }
  \item{at}{
    Character string specifying whether to compute the probability values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
  \item{weights}{
    Optional. Weights for the data points of \code{X}.
    \pointweightsVector.
  }
  \item{relative}{
    Logical.
    If \code{FALSE} (the default) the algorithm
    computes the probabilities of each type of point.
    If \code{TRUE}, it computes the    
    \emph{relative risk}, the ratio of probabilities
    of each type relative to the probability of a control.
  }
  \item{normalise}{
    Logical value specifying whether the results should be normalised
    so that constant risk corresponds to the value 1.
  }
  \item{casecontrol}{
    Logical. Whether to treat a bivariate point pattern
    as consisting of cases and controls, and return only the
    probability or relative risk of a case.
    Ignored if there are more than 2 types of points.
    See Details.
  }
  \item{control}{
    Integer, or character string, identifying which mark value
    corresponds to a control. 
  }
  \item{case}{
    Integer, or character string, identifying which mark value
    corresponds to a case (rather than a control)
    in a bivariate point pattern.
    This is an alternative to the argument \code{control}
    in a bivariate point pattern. 
    Ignored if there are more than 2 types of points.
  }
}
\details{
  The function \code{relriskAdaptiveKernel} is generic,
  with a method for point patterns of class \code{"ppp"}
  documented here.

  This function performs the same task as
  \code{\link[spatstat.explore]{relrisk.ppp}} except that
  \code{\link[spatstat.explore]{relrisk.ppp}} performs fixed-bandwidth kernel
  smoothing using a single value of smoothing bandwidth \code{sigma},
  whereas \code{relriskAdaptiveKernel.ppp} performs variable-bandwidth
  kernel smoothing using a different value of smoothing bandwidth for
  each point in the pattern \code{X}.

  Calculation of standard errors is not yet supported.
}
\value{
  If \code{X} consists of only two types of points,
  and if \code{casecontrol=TRUE},
  the result is a pixel image (if \code{at="pixels"})
  or a vector (if \code{at="points"}).
  The pixel values or vector values
  are the probabilities of a case if \code{relative=FALSE},
  or the relative risk of a case (probability of a case divided by the
  probability of a control) if \code{relative=TRUE}.

  If \code{X} consists of more than two types of points,
  or if \code{casecontrol=FALSE}, the result is:
  \itemize{
    \item (if \code{at="pixels"})
    a list of pixel images, with one image for each possible type of point.
    The result also belongs to the class \code{"solist"} so that it can
    be printed and plotted.
    \item
    (if \code{at="points"})
    a matrix of probabilities, with rows corresponding to
    data points \eqn{x_i}{x[i]}, and columns corresponding
    to types \eqn{j}.
  }
  The pixel values or matrix entries
  are the probabilities of each type of point if \code{relative=FALSE},
  or the relative risk of each type (probability of each type divided by the
  probability of a control) if \code{relative=TRUE}.

  If \code{relative=FALSE}, the resulting values always lie between 0
  and 1. If \code{relative=TRUE}, the results are either non-negative
  numbers, or the values \code{Inf} or \code{NA}.

  If \code{X} has \code{m > 1} columns of marks, the result is a list of
  length \code{m} containing the results for each column of marks
  as described above.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.explore]{relrisk.ppp}},
  \code{\link[spatstat.explore]{densityAdaptiveKernel.ppp}},
  \code{\link[spatstat.explore]{bw.abram.ppp}}.
}
\examples{
  bwa <- 2 * nndist(amacrine, proper=TRUE)
  A <- relriskAdaptiveKernel(amacrine, bw=bwa)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
