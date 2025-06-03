\name{auc}
\alias{auc}
\alias{auc.ppp}
\alias{auc.roc}
\alias{auc.cdftest}
\alias{auc.bermantest}
\alias{auc.im}
\title{
  Area Under ROC Curve for Point Pattern Data
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for a point pattern or other data.
}
\usage{
auc(X, \dots)

\method{auc}{ppp}(X, covariate, \dots, high = TRUE,  subset=NULL)

\method{auc}{roc}(X, \dots)

\method{auc}{cdftest}(X, \dots, high=TRUE)

\method{auc}{bermantest}(X, \dots, high=TRUE)

\method{auc}{im}(X, covariate, \dots, high=TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"lpp"})
    or fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"lppm"})
    or fitted spatial logistic regression model
    (object of class \code{"slrm"}) 
    or an ROC curve
    (object of class \code{"roc"} computed by \code{\link[spatstat.explore]{roc}}).
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.    
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{roc}},
    and arguments passed to \code{\link[spatstat.geom]{as.mask}}
    controlling the pixel resolution for calculations,
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, for which the AUC should be
    calculated.
  }
}
\details{
  This command computes the AUC, the area under the Receiver Operating
  Characteristic curve. The ROC itself is computed by \code{\link[spatstat.explore]{roc}}.

  The function \code{auc} is generic. There are methods for
  point patterns, fitted point process models, and many other kinds of
  objects.
  
  For a point pattern \code{X} and a covariate \code{Z}, the
  AUC is a numerical index that measures the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  Let \eqn{x_i}{x[i]} be a randomly-chosen data point from \code{X}
  and \eqn{U} a randomly-selected location in the study region.
  The AUC is the probability that
  \eqn{Z(x_i) > Z(U)}{Z(x[i]) > Z(U)}
  assuming \code{high=TRUE}.
  That is, AUC is the probability that a randomly-selected data point
  has a higher value of the covariate \code{Z} than does a
  randomly-selected spatial location. The AUC is a number between 0 and 1. 
  A value of 0.5 indicates a complete lack of discriminatory power.
  
  Methods for calculating AUC for a point process model or
  spatial logistic regression model are described in
  \code{\link[spatstat.model]{auc.ppm}} and \code{\link[spatstat.linnet]{auc.lpp}}.

  Some other kinds of objects in \pkg{spatstat} contain sufficient data to
  compute the AUC. These include the objects returned by
  \code{\link[spatstat.explore]{rhohat}},
  \code{\link[spatstat.explore]{cdf.test}} and \code{\link[spatstat.explore]{berman.test}}. Methods are
  provided here to compute the AUC from these objects.
}
\value{
  Numeric.
  For \code{auc.ppp}, \code{auc.cdftest}, \code{auc.bermantest}
  and \code{auc.im}, the result is a single number
  giving the AUC value.

  For \code{auc.roc}, the result is a numeric vector with one entry
  for each column of function values of \code{X}.
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
  \code{\link[spatstat.explore]{roc}},
  \code{\link[spatstat.model]{auc.ppm}},
  \code{\link[spatstat.linnet]{auc.lpp}}.

  \code{\link[spatstat.explore]{youden}}.
}
\examples{
  auc(swedishpines, "x")
}
\keyword{spatial}
