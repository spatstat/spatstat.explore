#'
#' tolcon.R
#' 
#' Generate tolerance contours for relative risk or conditional probability
#'
#' Copyright (c) 2019-2025 Adrian Baddeley and Tilman Davies
#' GNU Public Licence >= 2.0
#'
#' $Revision: 1.4 $ $Date: 2026/01/21 06:26:39 $

tolcon <- function(X, ..., nsim=19,
                   alternative=c("greater", "less", "two.sided"),
                   verbose=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  check.1.integer(nsim)
  alternative <- match.arg(alternative)
  Zdata <- as.solist(relrisk(X, ...))
  n <- length(Zdata)
  rankcount <- rep(list(0), n)
  if(verbose) 
    cat(paste("Simulating", nsim, "random labellings..."))
  pstate <- list()
  for(i in seq_len(nsim)) {
    Xsim <- rlabel(X)
    Zsim <- as.solist(relrisk(Xsim, ...))
    if(length(Zsim) != n)
      stop("Different numbers of images produced in data and simulation")
    for(j in seq_len(n)) 
      rankcount[[j]] <- rankcount[[j]] + (Zsim[[j]] >= Zdata[[j]])
    if(verbose)
      pstate <- progressreport(i, nsim, state=pstate)
  }
  if(verbose) 
    splat("Done.")
  result <- Zdata
  for(j in seq_len(n)) {
    y <- result[[j]]
    r <- rankcount[[j]]
    attr(y, "pvalues") <-
      switch(alternative,
             greater = { (r+1)/(nsim+1.0) },
             less    = { (nsim-r+1)/(nsim+1.0) },
             two.sided = {
               pg <- (r+1)/(nsim+1.0)
               pl <- (nsim-r+1)/(nsim+1.0)
               eval.im(pmin(1, 2*pmin(pg, pl)))
             })
    attr(y, "nsim") <- nsim
    attr(y, "alternative") <- alternative
    class(y) <- unique(c("tolcon", class(y)))
    result[[j]] <- y
  }
  class(result) <- unique(c("tolconlist", class(result)))
  return(result)
}

shift.tolcon <- function(X, ...) {
  y <- shift.im(X, ...)
  attr(y, "pvalues") <- shift.im(attr(X, "pvalues"), ...)
  return(y)
}

shift.tolconlist <- function(X, ...) {
  y <- solapply(X, shift, ...)
  class(y) <- union("tolconlist", class(y))  
  return(y)
}

plot.tolconlist <- function(x, ...) {
  plot.imlist(x, ..., plotcommand=plot.tolcon)
}

plot.tolcon <- function(x, ...,
                        show.contour=TRUE,
                        levels=0.05) {
  result <- plot.im(x, ...)
  if(show.contour)
    do.call.matched(contour.im,
                    resolve.defaults(list(x=attr(x, "pvalues"),
                                          levels=levels * 1.0001,
                                          add=TRUE),
                                     list(...),
                                     list(col="white",
                                          drawlabels=FALSE)),
                    extrargs=names(formals(contour.default)))
  return(result)
}

print.tolcon <- function(x, ...) {
  print.im(x, ...) 
  splat("\t[Includes tolerance p-values based on",
        attr(x, "nsim"),
        "random labellings]")
  return(invisible(NULL))
}
