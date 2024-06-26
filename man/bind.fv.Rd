\name{bind.fv}
\alias{bind.fv}
\alias{cbind.fv}
\title{
  Combine Function Value Tables
}
\description{
  Advanced Use Only.
  Combine objects of class \code{"fv"},
  or glue extra columns of data onto an existing \code{"fv"} object.
}
\usage{
\method{cbind}{fv}(...)

bind.fv(x, y, labl = NULL, desc = NULL, preferred = NULL, clip=FALSE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, which are objects of class \code{"fv"},
    or other data. See Details.
  }
  \item{x}{
    An object of class \code{"fv"}.
  }
  \item{y}{
    Either an object of class \code{"fv"}, a data frame,
    or a function. See Details.
  }
  \item{labl}{
    Plot labels (see \code{\link{fv}}) for columns of \code{y}.
    A character vector. 
  }
  \item{desc}{
    Descriptions (see \code{\link{fv}})
    for columns of \code{y}. A character vector.
  }
  \item{preferred}{
    Character string specifying the column which is to be the
    new recommended value of the function.
  }
  \item{clip}{
    Logical value indicating whether each object must have exactly the
    same domain, that is, the same sequence of values of the function argument
    (\code{clip=FALSE}, the default) or whether objects with different
    domains are permissible and will be restricted
    to a common domain (\code{clip=TRUE}).
  }  
}
\details{
  This documentation is provided
  for experienced programmers who want to modify the internal
  behaviour of \pkg{spatstat}.

  The function \code{cbind.fv} is a method for the generic
  \R function \code{\link{cbind}}. It combines any number of
  objects of class \code{"fv"} into a single object of
  class \code{"fv"}. The objects must be compatible, in the sense
  that they have identical values of the function argument.
  
  The function \code{bind.fv} is a lower level
  utility which glues additional columns onto an
  existing object \code{x} of class \code{"fv"}.
  It has three modes of use:
  \itemize{
    \item 
    If the additional dataset \code{y} is an object of class \code{"fv"}, then
    \code{x} and \code{y} must be compatible as described above.
    Then the columns of \code{y} that contain function values
    will be appended to the object \code{x}.
    \item 
    Alternatively if \code{y} is a data frame, then \code{y} must have the
    same number of rows as \code{x}. All columns of \code{y} will be
    appended to \code{x}.
    \item
    Alternatively if \code{y} is a function in the \R language, then
    this function will be evaluated at the argument values stored
    in the object \code{x}, and these function values will be appended
    as a new column to \code{x}.
  }
  The arguments \code{labl} and \code{desc} provide
  plot labels and description strings (as described in \code{\link{fv}})
  for the \emph{new} columns. If \code{y} is an object of class
  \code{"fv"} then \code{labl} and \code{desc} are optional, and
  default to the relevant entries in the object \code{y}. 
  If \code{y} is a data frame then
  \code{labl} and \code{desc} should be provided, but there is a
  default.

  For additional flexibility, \code{cbind.fv} also accepts arguments
  which are data frames or functions.
}
\value{
  An object of class \code{"fv"}.
}
\author{
  \spatstatAuthors.
}
\examples{
   K1 <- Kest(cells, correction="border")
   K2 <- Kest(cells, correction="iso")

   # remove column 'theo' to avoid duplication
   K2 <- K2[, names(K2) != "theo"]

   cbind(K1, K2)

   bind.fv(K1, K2, preferred="iso")

   # constrain border estimate to be monotonically increasing
   bm <- cumsum(c(0, pmax(0, diff(K1$border))))
   bind.fv(K1, data.frame(bmono=bm),
               "\%s[bmo](r)",
               "monotone border-corrected estimate of \%s",
               "bmono")

   # add a column of values defined by a function
   cbind(K1, upper=function(r) { pi * r^2 + 0.1 })
}
\seealso{
  \code{\link{fv}} for creating objects of class \code{"fv"} from raw data.

  \code{\link{collapse.fv}} for combining several \code{"fv"} objects
  with similar columns.
  
  \code{\link{with.fv}} for evaluating expressions.

  \code{\link{fvnames}} for extracting and assigning the column names
  of standard components of \code{"fv"} objects.

  \emph{Undocumented} functions for modifying an \code{"fv"} object
  include \code{tweak.fv.entry} and \code{rebadge.fv}.
}
\keyword{spatial}
\keyword{attribute}

