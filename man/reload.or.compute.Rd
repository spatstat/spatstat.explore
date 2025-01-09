\name{reload.or.compute}
\alias{reload.or.compute}
\title{
  Perform Computations or Retrieve Results From File
}
\description{
  This utility either performs computations and saves the results in
  a file, or retrieves the results of previous computations stored in a file.
  If the designated file does not yet exist,
  the expression will be evaluated, and the results will be saved in the file.
  If the file already exists, the results will be re-loaded from the file.
}
\usage{
reload.or.compute(filename, expr, objects = NULL,
                  context = parent.frame(), 
                  destination = parent.frame(),
                  force=FALSE, verbose=TRUE, exclude=NULL)
}
\arguments{
  \item{filename}{
    Name of data file. A character string.
  }
  \item{expr}{
    \R language expression to be evaluated.
  }
  \item{objects}{
    Optional character vector of names of objects to be saved
    in \code{filename} after evaluating \code{expr},
    or names of objects that should be present in \code{filename}
    when loaded.
  }
  \item{exclude}{
    Optional character vector of names of objects
    that should \bold{not} be saved in \code{filename}
    and are \bold{not} expected to be present in \code{filename}.
  }
  \item{context}{
    Environment containing objects that are mentioned in \code{expr}
    (other than objects in the global environment).
  }
  \item{destination}{
    Environment into which the resulting objects should be assigned.
  }
  \item{force}{
    Logical value indicating whether to perform the computation
    in any case.
  }
  \item{verbose}{
    Logical value indicating whether to print a message indicating
    whether the data were recomputed or reloaded from the file.
  }
}
\details{
  This facility is useful for saving, and later re-loading, the results of
  time-consuming computations. It would typically be
  used in an \R script file or an \code{\link[utils]{Sweave}} document.

  If the file called \code{filename} does not yet exist
  (or if \code{force=TRUE}),
  then \code{expr} will be evaluated
  and the results will be saved in \code{filename}
  using \code{\link[base]{save}}.
  By default, all objects that were created by evaluating the expression
  will be saved in the file.
  The optional argument \code{objects} specifies which results should be saved
  to the file. The optional argument \code{exclude}
  specifies results which should \emph{not} be saved.

  If the file called \code{filename} already exists
  (and if \code{force=FALSE}, the default), then this file
  will be loaded into \R using \code{\link{load}}.
  The optional argument \code{objects} specifies the names
  of objects that must be present in the file; a warning is issued
  if any of them are missing.

  The resulting objects (either evaluated or loaded from file)
  can be assigned into any desired
  \code{destination} environment.
  The default behaviour is equivalent to evaluating \code{expr}
  in the current environment.

  If \code{force=TRUE} then \code{expr} will be evaluated
  (regardless of whether the file already exists or not)
  and the results will be saved in \code{filename}, overwriting
  any previously-existing file with that name. This is a convenient
  way to force the code to re-compute everything
  in an \R script file or \code{\link[utils]{Sweave}} document.
}
\value{
  Character vector (invisible) giving the names of the objects computed
  or loaded.
}
\examples{
  ## Demonstration using a temporary file
  ## (For real applications, use a permanent file in your own filespace)
  myfile <- paste0(tempdir(), .Platform$file.sep, "mydata.rda")
  reload.or.compute(myfile, {
     # some very long computation ending with ..
     x <- 42
     intermediateWorking <- 12345
     y <- sqrt(x)
  }, exclude="intermediateWorking")
  ## the values x and y are saved
}
\author{
  \adrian
  and \rolf.
}
\keyword{utilities}
