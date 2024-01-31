#' boyce.R
#'
#' Discrete and continuous Boyce index
#'
#' $Revision: 1.3 $  $Date: 2024/01/31 06:59:18 $
#' 
#' Copyright (c) 2024 Adrian Baddeley

boyce <- function(X, Z, ..., breaks=NULL, halfwidth=NULL) {
  stopifnot(is.ppp(X))
  lbar <- intensity(unmark(X))
  if(is.im(Z) && Z$type == "factor") {
    ## convert to tessellation
    Z <- tess(image=Z)
  }
  if(is.tess(Z)) {
    ## discrete Boyce index
    Y <- as.tess(quadratcount(X, tess=Z))
    lam <- marks(Y)[,1]/tile.areas(Y)
    result <- Z
    marks(result) <- lam/lbar
  } else {
    ## continuous Boyce index
    ngiven <- (!is.null(breaks)) + (!is.null(halfwidth))
    if(ngiven == 0) stop("Either 'breaks' or 'halfwidth' should be given")
    if(ngiven == 2) stop("Arguments 'breaks' and 'halfwidth' are incompatible")
    if(!is.null(breaks)) {
      ## discrete Boyce index based on bands of Z values
      result0 <- rhohat(X, Z, smoother="piecewise", breaks=breaks,
                        ...)
    } else {
      ## continuous Boyce index, equivalent to rhohat using box kernel
      result0 <- rhohat(X, Z, kernel="rect", bw=halfwidth/kernel.factor("rect"),
                        ...)
    }
    result <- eval.fv(result0/lbar)
    newylab <- substitute(B(x), list(x=as.name(fvnames(result0, ".x"))))
    oldlabels <- attr(result0, "labl")
    newlabels <- oldlabels[colnames(result0) %in% colnames(result)]
    result <- rebadge.fv(result, new.fname="B", new.ylab=newylab,
                         new.labl=newlabels)
  }
  return(result)
}

