\name{roc}
\alias{roc}
\alias{roc.ppp}
\alias{roc.cdftest}
\alias{roc.bermantest}
\alias{roc.im}
\title{
  Receiver Operating Characteristic
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern or a fitted point process model.
}
\usage{
 roc(X, \dots)

\method{roc}{ppp}(X, covariate, 
                  \dots,
                  baseline = NULL, high = TRUE, weights = NULL,
                  observations=c("exact", "presence"),
                  method = "raw",
                  CI = "none", alpha=0.05,
                  subset=NULL)

\method{roc}{cdftest}(X, \dots, high=TRUE)

\method{roc}{bermantest}(X, \dots, high=TRUE)

\method{roc}{im}(X, covariate, \dots, high=TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"lpp"})
    or fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"lppm"})
    or fitted spatial logistic regression model
    (object of class \code{"slrm"}) or some other kind of data.
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.
    Traditionally omitted when \code{X} is a fitted model.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} controlling the
    pixel resolution for calculations.
  }
  \item{baseline}{
    Optional. A spatial object giving a baseline intensity.
    Usually a \code{function(x,y)} or
    a pixel image (object of class \code{"im"})
    giving the baseline intensity at
    any location within the observation window.
    Alternatively a point pattern (object of class \code{"ppp"})
    with the locations of the reference population.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
  \item{weights}{
    Optional. Numeric vector of weights attached to the data points.
  }
  \item{observations}{Character string (partially matched)
    specifying whether to compute the ROC curve using the
    exact point coordinates (\code{observations="exact"}, the default)
    or using the discretised presence-absence data
    (\code{observations="presence"}).
  }
  \item{method}{
    The method or methods that should be used to estimate the ROC curve.
    A character vector: current choices are
    \code{"raw"}, \code{"monotonic"}, \code{"smooth"} and \code{"all"}.
    See Details.
  }
  \item{CI}{
    Character string (partially matched) specifying whether confidence
    intervals should be computed, and for which method.
    See Details.
  }
  \item{alpha}{
    Numeric value between 0 and 1. The confidence intervals will have
    confidence level \code{1-alpha}. The default gives 95\%
    confidence intervals.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, from which the ROC should be
    calculated.
  }
}
\details{
  This command computes the Receiver Operating Characteristic (ROC)
  curve. The area under the ROC is computed by \code{\link{auc}}.

  The function \code{roc} is generic, with methods for point patterns,
  fitted point process models, and other kinds of data.
  
  For a point pattern \code{X} and a spatial covariate \code{Z}, the
  ROC is a plot showing the ability of the
  \code{covariate} to separate the spatial domain
  into areas of high and low density of points.
  For each possible threshold \eqn{z}, the algorithm calculates
  the fraction \eqn{a(z)} of area in the study region where the
  covariate takes a value greater than \eqn{z}, and the
  fraction \eqn{b(z)} of data points for which the covariate value
  is greater than \eqn{z}. The ROC is a plot of \eqn{b(z)} against
  \eqn{a(z)} for all thresholds \eqn{z}. This is called the \sQuote{raw}
  ROC curve.

  There are currently three methods to estimate the ROC curve:
  \describe{
    \item{\code{"raw"}}{
      uses the raw empirical spatial cummulative distribution function of the
      covariate.
    }
    \item{\code{"monotonic"}}{
      uses a monotonic regression to estimate the relation between the covariate
      and the point process intensity and then calculates the ROC from that.
      This corresponds to a either a convex minorant or a concave majorant of
      the raw ROC curve.
    }
    \item{\code{"smooth"}}{
      uses a smooth estimate of the relation between the covariate and the point
      process intensity and then calculates the ROC from that. See
      \code{\link{roc.rhohat}} for details.
    }
    \item{\code{"all"}}{
      uses all of the above methods.
    }
  }

  If \code{CI} is one of the strings \code{'raw'},
  \code{'monotonic'} or \code{'smooth'}, then
  pointwise 95\% confidence intervals for the true ROC curve
  will be computed based on the \code{raw}, \code{monotonic} or
  \code{smooth} estimates, respectively.
  The confidence level is \code{1-alpha}, so that for example
  \code{alpha=0.01} would give 99\% confidence intervals.
  By default, confidence bands for the ROC curve are not computed.
  
  Some other kinds of objects in \pkg{spatstat} contain sufficient data to
  compute the ROC curve. These include the objects returned by
  \code{\link[spatstat.explore]{rhohat}},
  \code{\link[spatstat.explore]{cdf.test}} and \code{\link[spatstat.explore]{berman.test}}. Methods are
  provided here to compute the ROC curve from these objects.

  The method for pixel images (objects of class \code{"im"})
  assumes that \code{X} represents a density or intensity function,
  and that the objective is to segregate the spatial region into
  subregions of high and low total density by thresholding the
  \code{covariate}.
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
  Also belongs to class \code{"roc"}.
}
\references{
  \rocpaper.

  Lobo, J.M.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S.
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \ifelse{latex}{\out{Birkh{\"a}user}}{Birkhauser}, Basel.
}
\author{
  \rocketAuthors.
}
\seealso{
  \code{\link[spatstat.model]{roc.ppm}},
  \code{\link[spatstat.linnet]{roc.lpp}},
  \code{\link[spatstat.explore]{roc.rhohat}}.

  \code{\link[spatstat.explore]{auc}}
}
\examples{

  gold <- rescale(murchison$gold, 1000, "km")
  faults <- rescale(murchison$faults, 1000, "km")
  dfault <- distfun(faults)

  if(interactive()) {
    plot(roc(gold, dfault, method = "all", high=FALSE))
  } else {
    ## reduce sample resolution to save computation time in test
    plot(roc(gold, dfault, method = "all", high=FALSE, eps=8))
  }

  # Using either an image or reference population as baseline
  cases <- split(chorley)$larynx
  controls <- split(chorley)$lung
  covar <- distfun(as.ppp(chorley.extra$incin, W = Window(chorley)))
  if(interactive()) {
    population <- density(controls, sigma=0.15, eps=0.1)
  } else {
    ## reduce resolution to save computation time in test
    population <- density(controls, sigma=0.3, eps=0.25)
  }
  population <- eval.im(pmax(population, 1e-10))
  roc1 <- roc(cases, covar, baseline = population, high = FALSE, method="all")
  roc2 <- roc(cases, covar, baseline = controls, high = FALSE, method="all")
  plot(anylist(roc1=roc1, roc2=roc2), main = "")
}
\keyword{spatial}
\keyword{nonparametric}
