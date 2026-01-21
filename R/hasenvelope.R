#'
#'    hasenvelope.R
#'
#'    A simple class of objects which contain additional envelope data
#' 
#'    $Revision: 1.2 $ $Date: 2026/01/21 06:26:39 $

hasenvelope <- function(X, E=NULL) {
  if(inherits(E, "envelope")) {
    attr(X, "envelope") <- E
    class(X) <- unique(c("hasenvelope", class(X)))
  }
  return(X)
}

print.hasenvelope <- function(x, ...) {
  NextMethod("print")
  splat("[Object contains simulation envelope data]")
  return(invisible(NULL))
}

envelope.hasenvelope <- function(Y, ..., Yname=NULL) {
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  E <- attr(Y, "envelope")
  return(envelope(E, ..., Yname=Yname))
}
  
