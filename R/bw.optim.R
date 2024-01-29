#
# bw.optim.R
#
#  Class of optimised bandwidths
#  Plotting the object displays the optimisation criterion
#
#  $Revision: 1.37 $  $Date: 2024/01/29 07:09:03 $
#

bw.optim <- function(cv, h,
                     iopt=if(optimum == "min") which.min(cv) else which.max(cv),
                     ...,
                     cvname, hname,
                     criterion="cross-validation",
                     optimum = c("min", "max"), 
                     warnextreme=TRUE, hargnames=NULL,
                     yexp=NULL,
                     unitname=NULL,
                     template=NULL,
                     exponent=1,
                     hword) {
  if(missing(cvname) || is.null(cvname)) cvname <- short.deparse(substitute(cv))
  if(missing(hname) || is.null(hname)) hname <- short.deparse(substitute(h))
  stopifnot(is.numeric(cv))
  stopifnot(is.numeric(h))
  stopifnot(length(h) == length(cv))
  optimum <- match.arg(optimum)
  result <- h[iopt]
  if(warnextreme) {
    optimised <- switch(optimum, min="minimised", max="maximised")
    if(is.infinite(result)) {
      warning(paste(criterion, "criterion was", optimised, "at",
                    hname, "=", as.numeric(result)),
              call.=FALSE)
    } else if((iopt == length(h) || iopt == 1)) {
      warning(paste(criterion, "criterion was", optimised, "at",
                    if(iopt == 1) "left-hand" else "right-hand",
                    "end of interval",
                    paste0(prange(signif(range(h[is.finite(h)]), 3)), ";"), 
                    "use", ngettext(length(hargnames), "argument", "arguments"),
                    paste(sQuote(hargnames), collapse=", "),
                    "to specify a wider interval for bandwidth", sQuote(hname)),
              call.=FALSE)
    }
  }
  if(missing(hword))
    hword <- if(is.null(template)) "bandwidth" else "scale factor"
  attr(result, "cv") <- cv
  attr(result, "h") <- h
  attr(result, "iopt") <- iopt
  attr(result, "labels") <- list(hname=hname, cvname=cvname)
  attr(result, "info") <- list(...)
  attr(result, "criterion") <- criterion
  attr(result, "optimum") <- optimum
  attr(result, "hargnames") <- hargnames
  attr(result, "units") <- as.unitname(unitname)
  attr(result, "yexp") <- yexp
  attr(result, "template") <- template
  attr(result, "exponent") <- exponent %orifnull% 1
  attr(result, "hword") <- hword
  class(result) <- "bw.optim"
  return(result)
}

print.bw.optim <- function(x, ...) {
  y <- as.numeric(x)
  names(y) <- hname <- attr(x, "labels")$hname
  print(y, ...)
  if(!is.null(m <- attr(x, "template"))) {
    exponent <- attr(x, "exponent") %orifnull% 1
    hpow <- if(exponent == 1) hname else paste0(paren(hname), "^", exponent)
    cat("\n")
    splat(hpow, "is interpreted as a multiple of:")
    print(m)
  }
  return(invisible(NULL))
}

summary.bw.optim <- function(object, ...) {
  z <- attributes(object)
  z$hopt <- hopt <- as.numeric(object)
  z$is.extreme <- is.infinite(hopt) || with(z, iopt == 1 || iopt == length(h))
  structure(z, class="summary.bw.optim")
}

print.summary.bw.optim <- function(x, ..., digits=3) {
  hword <- x$hword %orifnull% "bandwidth"
  Hword <- paste0(toupper(substring(hword, 1, 1)), substring(hword, 2))
  splat(Hword, "value selected by", x$criterion)
  su <- summary(x$units)
  splat("Optimal value:",
        x$labels$hname, "=",
        signif(x$hopt, digits=digits),
        if(x$hopt == 1) su$singular else su$plural,
        su$explain)
  splat("Search performed over", length(x$h),
        "candidate values of", hword, "in the interval",
        prange(signif(range(x$h), digits=digits)))
  optname <- if(is.null(x$optimum)) "Optimum" else 
             switch(x$optimum, min="Minimum", max="Maximum", x$optimum)
  splat(optname, "value of criterion",
        paste0(x$labels$cvname, ":"),
        signif(x$cv[x$iopt], digits=digits))
  if(isTRUE(x$is.extreme)) {
    splat(optname, "achieved at",
          if(is.infinite(x$hopt)) "infinity" else
          if(x$iopt == 1) "lower limit of range" else "upper limit of range")
  }
  if(!is.null(creator <- x$info$creator)) 
    splat("Computed by the function", sQuote(creator))
  if(!is.null(tem <- x$template)) {
    exponent <- x$exponent
    Hpow <- if(exponent == 1) Hword else paste0(paren(Hword), "^", exponent)
    splat(Hpow, "is interpreted as a multiplier of:")
    print(tem)
  }
  return(invisible(NULL))
}

as.data.frame.bw.optim <- function(x, ...) {
  h <- attr(x, "h")
  cv <- attr(x, "cv")
  df <- data.frame(h, cv)
  labels <- attr(x, "labels")
  colnames(df) <- labels[c("hname", "cvname")]
  info <- attr(x, "info")
  if(length(info) > 0) {
    lenfs <- lengths(info)
    if(any(ok <- (lenfs == nrow(df)))) {
      df <- cbind(df, as.data.frame(info[ok]))
    }
  }
  return(df)
}

as.fv.bw.optim <- function(x) {
  # convert to fv object
  df <- as.data.frame(x)
  dfnames <- colnames(df)
  hname <- dfnames[1L]
  cvname <- dfnames[2L]
  descrip <- c("smoothing parameter",
               paste(attr(x, "criterion"), "criterion"))
  if(ncol(df) > 2)
    descrip <- c(descrip, paste("Additional variable", sQuote(dfnames[-(1:2)])))
  labl <- c(hname, paste0(dfnames[-1L], paren(hname)))
  yexp <- attr(x, "yexp") %orifnull% substitute(CV(h),
                                                list(CV=as.name(cvname),
                                                     h=as.name(hname)))
  xfv <- fv(df,
            argu=hname,
            ylab=yexp,
            valu=cvname,
            labl=labl,
            desc=descrip,
            fname=cvname,
            yexp=yexp)
  fvnames(xfv, ".") <- cvname
  unitname(xfv) <- unitname(x)
  return(xfv)
}

plot.bw.optim <- function(x, ...,
                          showopt=TRUE, optargs=list(lty=3, col="blue")) {
  xname <- short.deparse(substitute(x))
  # convert to fv object
  xfv <- as.fv(x)
  # plot cross-validation criterion
  out <- do.call(plot.fv,
                 resolve.defaults(list(x=quote(xfv)),
                                  list(...),
                                  list(main=xname)))
  # Turn off 'showopt' if the x-variable is not the bandwidth
  if(missing(showopt)) {
    argh <- list(...)
    isfmla <- unlist(lapply(argh, inherits, what="formula"))
    if(any(isfmla)) {
      fmla <- argh[[min(which(isfmla))]]
      xvar <- deparse(rhs.of.formula(fmla, tilde=FALSE))
      if(!(identical(xvar, fvnames(xfv, ".x")) || identical(xvar, ".x")))
        showopt <- FALSE
    }
  }
  # show optimal value?
  if(showopt) {
    hoptim <- as.numeric(x)
    if(spatstat.options('monochrome'))
      optargs <- col.args.to.grey(optargs)
    do.call(abline, append(list(v=hoptim), optargs))
  }
  if(is.null(out)) return(invisible(NULL))
  return(out)
}


