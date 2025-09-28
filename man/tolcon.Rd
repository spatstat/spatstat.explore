\name{tolcon}
\alias{tolcon}
\title{
  Tolerance Contours for Relative Risk
}
\description{
  Compute tolerance contours for the estimate of
  spatially-varying relative risk
  or spatially-varying probabilities of types.
}
\usage{
tolcon(X, \dots, nsim = 19, alternative=c("greater", "less", "two.sided"),
       verbose = TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern
    (object of class \code{"ppp"} with
    marks which are a factor).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{relrisk.ppp}}.
  }
  \item{nsim}{
    Number of simulations (random re-labellings of \code{X})
    to be used to compute the tolerance contours.
  }
  \item{alternative}{
    Character string (partially matched) specifying the direction in
    which values are deemed significant.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  This function computes tolerance contours for the
  images computed by \code{\link[spatstat.explore]{relrisk}}.
  
  The function \code{\link[spatstat.explore]{relrisk}} computes estimates of
  spatially-varying relative risk in a case-control study,
  or spatially-varying type distribution in a multitype point pattern.
  Tolerance contours for relative risk (Hazelton and Davies, 2009)
  are curves drawn around the regions where the estimated relative risk
  (or estimated probability of a given type) is significantly different from the
  average risk or average proportion.

  Significance is assessed by a Monte Carlo test. First the
  original dataset is analysed by calling
  \code{\link[spatstat.explore]{relrisk}(X, \dots)}. Then
  \code{X} is randomly relabelled (that is, the marks attached to the
  points are randomly permuted using \code{\link[spatstat.random]{rlabel}}) and the
  estimate of relative risk or type probability is computed from this
  relabelled data. This is repeated \code{nsim} times, yielding
  \code{nsim} relative risk images or \code{nsim} lists of images
  of the probabilities of each time. A Monte Carlo
  p-value is computed at each pixel; these values are multiples of
  \code{1/(nsim+1)}, so the default \code{nsim=19} produces p-values
  which are multiples of \code{1/20 = 0.05}.

  The result returned by \code{tolcon} is identical to the
  result returned by \code{\link[spatstat.explore]{relrisk}}, except that each of the
  images of relative risk or type probability belongs to the
  additional class \code{"tolcon"}. Each image of class \code{"tolcon"}
  has an attribute \code{"pvalues"} which is an image
  containing the relevant p-values. The class \code{"tolcon"} has a
  plot method which first plots the image and then draws the
  tolerance contour.
}
\value{
  An image or a list of images.
  Each image belongs to the additional class \code{"tolcon"}
  (for which there is a plot method) and
  has an attribute \code{"pvalues"} which is another image
  containing the p-values.
}
\references{
  Hazelton, M.L. and Davies, T.M. (2009)
  Inference based on kernel estimates of the relative
  risk function in geographical epidemiology.
  \emph{Biometrical Journal} \bold{51}, 98--109.
}
\author{
  \adrian, \tilman and \Lucia Cobo \Sanchez.
}
\seealso{
  \code{\link[spatstat.explore]{relrisk.ppp}}
}
\examples{
   Z <- tolcon(mucosa, casecontrol=FALSE, sigma=0.1)
   plot(Z, main="")
}
\keyword{spatial}
\keyword{smooth}
\keyword{htest}
