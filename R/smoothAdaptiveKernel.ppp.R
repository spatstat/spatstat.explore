#'
#'   smoothAdaptiveKernel.ppp.R
#'
#'   $Revision: 1.2 $  $Date: 2026/07/15 04:45:51 $
#'
#'
#'  Adaptive kernel smoothing via 3D FFT
#'

SmoothAdaptiveKernel <- function(X, bw, ...) {
  UseMethod("SmoothAdaptiveKernel")
}

SmoothAdaptiveKernel.ppp <- function(X, bw, ...,
                                     weights=NULL,
                                     at=c("pixels", "points")) {
  stopifnot(is.ppp(X))
  at <- match.arg(at)
  if(!is.marked(X, dfok=TRUE, na.action="fatal"))
    stop("X should be a marked point pattern", call.=FALSE)
  X <- coerce.marks.numeric(X)
  marx <- marks(X)
  if(!all(is.finite(as.matrix(marx))))
    stop("Some mark values are Inf, NaN or NA", call.=FALSE)
  univariate <- is.null(dim(marx))
  nc <- if(univariate) 1 else ncol(marx)
  nX <- npoints(X)

  ## trivial case
  if(nX == 0) {
    cn <- colnames(marks(X))
    switch(at,
           points = {
             result <- if(univariate) numeric(0) else
                       matrix(, 0, nc, dimnames=list(NULL, cn))
           },
           pixels = {
             result <- as.im(NA_real_, Window(X))
             if(!univariate) {
               result <- as.solist(rep(list(result), nc))
               names(result) <- cn
             }
           })
    return(result)
  }

  ## ensure weights are numeric
  if(weighted <- !missing(weights) && !is.null(weights)) {
    pa <- parent.frame()
    weights <- pointweights(X, weights=weights, parent=pa)
    weighted <- !is.null(weights)
  } else weights <- NULL

  ## calculate denominator and determine bandwidths
  denom <- densityAdaptiveKernel(unmark(X), bw=bw, ..., at=at,
                                 weights=weights)
  bw <- attr(denom, "bw")

  ## prepare nearest neighbour information
  nn <- switch(at,
               points = nnwhich(X),
               pixels = nnmap(X, what="which"))

  ## calculate numerator(s)
  results <- vector(mode="list", length=nc)
  if(univariate) {
    weightsmarks <- if(weighted) (weights * marx) else marx
    numer <- densityAdaptiveKernel(unmark(X), bw=bw, ..., at=at,
                                   weights=weightsmarks)
    results[[1L]] <- saferatio(numer, denom, marx, nn)
  } else {
    for(j in 1:nc) {
      weightsmarksj <- if(weighted) (weights * marx[,j]) else marx[,j]
      numerj <- densityAdaptiveKernel(unmark(X), bw=bw, ..., at=at,
                                   weights=weightsmarksj)
      results[[j]] <- saferatio(numerj, denom, marx, nn)
    }
    names(results) <- cn
  }
  ## format of result
  result <- if(univariate) results[[1L]] else
            switch(at, 
                   points = do.call(cbind, results),
                   pixels = as.solist(results))
  attr(result, "bw") <- bw
  return(result)
}

saferatio <- function(numer, denom, values, nn) {
  eps <- .Machine$double.eps
  if(is.numeric(numer) && is.numeric(denom)) {
    ratio <- numer/denom
    ## trap small values of denominator
    nbg <- !is.finite(ratio) | (denom < eps)
    if(any(nbg)) {
      warning(paste("Very small denominator detected:",
                    "sigma is probably too small"),
              call.=FALSE)
      ## l'Hopital's rule
      ratio[nbg] <- values[nn[nbg]]
      ## save warnings
      uhoh <- attr(numer, "warnings")
      attr(ratio, "warnings") <- unique(c(uhoh, "underflow"))
    }
  } else if(is.im(numer) && is.im(denom)) {
    ratio <- eval.im(numer/denom)
    ## trap small values of denominator
    nbg <- eval.im(is.infinite(ratio)
                   | is.nan(ratio)
                   | (denom < eps))
    if(bad <- any(as.matrix(nbg), na.rm=TRUE)) {
      warning(paste("Numerical underflow detected:",
                    "sigma is probably too small"),
              call.=FALSE)
      ## l'Hopital's rule
      ratio[nbg] <- values[nn[nbg]]
      ## save warnings
      uhoh <- attr(numer, "warnings")
      attr(ratio, "warnings") <- unique(c(uhoh, "underflow"))
    }
  } else stop("Internal error: unrecognised format of ratio", call.=FALSE)
  return(ratio)
}
