##
## densityfun.R
##
## Exact 'funxy' counterpart of density.ppp
##
##  $Revision: 1.15 $ $Date: 2023/04/01 02:25:53 $


densityfun <- function(X, ...) {
  UseMethod("densityfun")
}

densityfun.ppp <- function(X, sigma=NULL, ...,
                          weights=NULL, edge=TRUE, diggle=FALSE) {
  verifyclass(X, "ppp")
  ## standard errors are not yet supported in densitycrossEngine
  if(isTRUE(list(...)$se))
    warning("Standard errors are not yet supported in densityfun.ppp",
            call.=FALSE)
  ## handle weights now
  weightsgiven <- !missing(weights) && !is.null(weights) 
  if(weightsgiven) {
    # convert to numeric
    if(is.im(weights)) {
      weights <- safelookup(weights, X) # includes warning if NA
    } else if(is.expression(weights)) 
      weights <- eval(weights, envir=as.data.frame(X), enclos=parent.frame())
    if(length(weights) == 0)
      weightsgiven <- FALSE
  }
  if(weightsgiven) {
    check.nvector(weights, npoints(X), vname="weights")
  } else weights <- NULL
  ## 
  stuff <- list(Xdata=X, weights=weights, edge=edge, diggle=diggle, ...)
  ##
  ## determine smoothing parameters
  ker <- resolve.2D.kernel(sigma=sigma, ...,
                           x=X, bwfun=bw.diggle, allow.zero=TRUE)
  stuff[c("sigma", "varcov")]  <- ker[c("sigma", "varcov")]
  ##
  g <- function(x, y=NULL, drop=TRUE) {
    Y <- xy.coords(x, y)[c("x", "y")]
    W <- Window(stuff$Xdata)
    ok <- inside.owin(Y, w=W)
    allgood <- all(ok)
    if(!allgood) Y <- lapply(Y, "[", i=ok)
    Xquery <- as.ppp(Y, W)
    vals <- do.call(densitycrossEngine, append(list(Xquery=Xquery), stuff))
    if(allgood || drop) return(vals)
    ans <- numeric(length(ok))
    ans[ok] <- vals
    ans[!ok] <- NA
    attr(ans, "sigma") <- attr(vals, "sigma")
    return(ans)
  }
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("densityfun", class(g))
  return(g)
}

print.densityfun <- function(x, ...) {
  cat("function(x,y)", "which returns",
      "kernel estimate of intensity for", fill=TRUE)
  X <- get("X", envir=environment(x))
  print(X, ...)
  cat("Optional argument:", "drop=TRUE", fill=TRUE)
  return(invisible(NULL))
}

## Method for as.im
## (enables plot.funxy, persp.funxy, contour.funxy to work for this class)

as.im.densityfun <- function(X, W=Window(X), ..., approx=TRUE) {
  if(!approx) {
    #' evaluate exactly at grid points using as.im.funxy -> as.im.function
    result <- as.im.function(X, W=W, ...)
  } else {
    #' faster, approximate evaluation using FFT
    stuff <- get("stuff", envir=environment(X))
    Xdata <- stuff[["Xdata"]]
    otherstuff <- stuff[names(stuff) != "Xdata"]
    if(!missing(W)) Xdata <- Xdata[W]
    result <- do.call(density,
                      resolve.defaults(list(x=quote(Xdata)),
                                       list(...), otherstuff))
  }
  return(result)
}


  
