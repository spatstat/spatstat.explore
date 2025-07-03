\name{cov.im}
\alias{cov.im}
\alias{cor.im}
\title{
  Covariance and Correlation between Images
}
\description{
  Compute the covariance or correlation between (the corresponding pixel
  values in) several images.
}
\usage{
cov.im(\dots, use = "complete.obs", method = c("pearson", "kendall", "spearman"))
cor.im(\dots, use = "complete.obs", method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is 
    a pixel image (object of class \code{"im"}).
    Alternatively, a single argument which is a list of pixel images.
  }
  \item{use}{
    Argument passed to \code{\link[stats]{cov}}
    or \code{\link[stats]{cor}} determining how to handle
    \code{NA} values in the data.
  }
  \item{method}{
    Argument passed to \code{\link[stats]{cov}}
    or \code{\link[stats]{cor}} determining the type of correlation
    that will be computed.
  }
}
\details{
  The arguments \code{\dots}
  should be pixel images (objects of class \code{"im"}).
  Their spatial domains must overlap, but need not have the same
  pixel dimensions.

  These functions compute the covariance or correlation between
  the corresponding pixel values in the images given.

  The pixel images are converted to a common pixel resolution (by resampling).
  Then the corresponding pixel values of each image are
  extracted. Finally the correlation or covariance between the
  pixel values of each pair of images, at corresponding pixels, is
  computed.

  The result is a symmetric matrix with one row and column for each image.
  The \code{[i,j]} entry is the correlation or covariance between
  the \code{i}th and \code{j}th images in the argument list.
  The row names and column names of the matrix are copied from the
  argument names if they were given (i.e. if the arguments were given
  as \code{name=value}).

  The argument \code{use} specifies how to handle \code{NA} values.
  A pixel value of \code{NA} is assigned to any
  pixel falling outside the spatial domain of an image
  (i.e. outside the window in which the image is defined).
  If any one of the image arguments \code{\dots} is defined on a
  non-rectangular window, or if the image arguments are not all defined
  on the same window, then the data will contain \code{NA} values.
  Options for the argument \code{use} are documented in the help file
  for \code{\link[stats]{cov}} and \code{\link[stats]{cor}}.
  \itemize{
    \item \code{use="complete.obs"} (the default):
    calculations are based on those pixels which lie inside the intersection
    of the windows of all the images. An error occurs if the
    intersection is empty.
    \item \code{use="na.or.complete"}:
    calculations are based on those pixels which lie inside the intersection
    of the windows of all the images. If the
    intersection is empty, a matrix of \code{NA} values is returned.
    \item \code{use="pairwise.complete.obs"}: the calculation of the
    covariance or correlation between each pair of images is
    based on the pixels which lie in the intersection of the
    windows of those two images. Only available when \code{method="pearson"}.
    The resulting matrix may not be positive definite.
    \item \code{use="everything"}: the calculation is based on all
    pixels, but any calculation of variance or covariance or correlation
    that includes an \code{NA} value gives an \code{NA} result in the
    corresponding entry in the matrix.
    \item \code{use="all.obs"}: the calculation is based on all pixels,
    and an error occurs if any pixel has an \code{NA} value in any image.
  }
  
  Note that \code{\link[stats]{cor}} and \code{\link[stats]{cov}}
  are not generic, so you have to type \code{cor.im}, \code{cov.im}.
}
\value{
  A symmetric matrix.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[stats]{cor}},
  \code{\link[stats]{cov}}

  \code{\link{pairs.im}}
}
\examples{
  cor.im(bei.extra)
}
\keyword{spatial}
\keyword{univar}
\keyword{nonparametric}

