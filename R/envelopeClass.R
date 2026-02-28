#'
#' envelopeClass.R
#'
#' support for class 'envelope'
#'
#' Methods for plot, print, summary, envelope, pool, as.data.frame
#' 
#'   $Revision: 1.2 $  $Date: 2026/02/28 12:22:02 $


as.data.frame.envelope <- function(x, ..., simfuns=FALSE) {
  if(simfuns && !is.null(sf <- attr(x, "simfuns"))) {
    # tack on the simulated functions as well
    y <- as.data.frame(bind.fv(x, sf, clip=TRUE))
    return(y)
  } 
  NextMethod("as.data.frame")
}

plot.envelope <- function(x, ..., main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  shade.given <- ("shade" %in% names(list(...)))
  shade.implied <- !is.null(fvnames(x, ".s"))
  if(!(shade.given || shade.implied)) {
    # ensure x has default 'shade' attribute
    # (in case x was produced by an older version of spatstat)
    if(all(c("lo", "hi") %in% colnames(x)))
      fvnames(x, ".s") <- c("lo", "hi")
    else warning("Unable to determine shading for envelope")
  }
  NextMethod("plot", main=main)
}

print.envelope <- function(x, ...) {
  e <- attr(x, "einfo")
  g <- e$global
  csr <- e$csr
  pois <- e$pois
  if(is.null(pois)) pois <- csr
  simtype <- e$simtype
  constraints <- e$constraints
  nr <- e$nrank
  nsim <- e$nsim
  V <- e$VARIANCE
  fname <- flat.deparse(attr(x, "ylab"))
  type <- if(V) paste("Pointwise", e$nSD, "sigma") else
          if(g) "Simultaneous" else "Pointwise"
  splat(type, "critical envelopes for", fname,
        "\nand observed value for", sQuote(e$Yname))
  if(!is.null(valname <- e$valname) && waxlyrical('extras'))
    splat("Edge correction:", dQuote(valname))
  ## determine *actual* type of simulation
  descrip <- e$sim.realisations %orifnull% 
    if(csr) "simulations of CSR"
    else if(!is.null(simtype)) {
      switch(simtype,
             csr="simulations of CSR",
             rmh=paste("simulations of fitted",
               if(pois) "Poisson" else "Gibbs",
               "model"),
             kppm="simulations of fitted cluster model",
             slrm="simulations of fitted spatial logistic regression model",
             expr="evaluations of user-supplied expression",
             func="evaluations of user-supplied function",
             list="point pattern datasets in user-supplied list",
             funs="columns of user-supplied data")
    } else "simulations of fitted model"
  if(!is.null(constraints) && nzchar(constraints))
    descrip <- paste(descrip, constraints)
  #  
  splat("Obtained from", nsim, descrip)
  #
  if(waxlyrical('extras')) {
    dual <- isTRUE(e$dual)
    usetheory <- isTRUE(e$use.theory)
    hownull <- if(usetheory) {
                 "(known exactly)"
               } else if(dual) {
                 paste("(estimated from a separate set of",
                       e$nsim2, "simulations)")
               } else NULL
    formodel <- if(csr) "for CSR" else NULL
    if(g) {
      splat("Envelope based on maximum deviation of", fname,
            "from null value", formodel, hownull)
    } else if(dual) {
      splat("Null value of", fname, formodel, hownull)
    }
    if(!is.null(attr(x, "simfuns"))) 
      splat("(All simulated function values are stored)")
    if(!is.null(attr(x, "simpatterns"))) 
      splat("(All simulated point patterns are stored)")
  }
  splat("Alternative:", e$alternative)
  if(!V && waxlyrical('extras')) {
    ## significance interpretation!
    alpha <- if(g) { nr/(nsim+1) } else { 2 * nr/(nsim+1) }
    splat("Significance level of",
          if(g) "simultaneous" else "pointwise",
          "Monte Carlo test:",
          paste0(if(g) nr else 2 * nr,
                 "/", nsim+1),
          "=", signif(alpha, 3))
  }
  if(waxlyrical('gory') && !is.null(pwrong <- attr(x, "pwrong"))) {
    splat("\t[Estimated significance level of pointwise excursions:",
          paste0("pwrong=", signif(pwrong, 3), "]"))
  }
  NextMethod("print")
}
                  
summary.envelope <- function(object, ...) {
  e <- attr(object, "einfo")
  g <- e$global
  V <- e$VARIANCE
  nr <- e$nrank
  nsim <- e$nsim
  csr <- e$csr
  pois <- e$pois
  if(is.null(pois)) pois <- csr
  simtype <- e$simtype
  constraints <- e$constraints
  alternative <- e$alternative
  use.theory <- e$use.theory
  has.theo <- "theo" %in% fvnames(object, "*")
  csr.theo <- csr && has.theo
  use.theory <- if(is.null(use.theory)) csr.theo else (use.theory && has.theo)
  fname <- short.deparse(attr(object, "ylab"))
  type <- if(V) paste("Pointwise", e$nSD, "sigma") else
          if(g) "Simultaneous" else "Pointwise"
  splat(type, "critical envelopes for", fname, 
      "\nand observed value for", sQuote(e$Yname))
  # determine *actual* type of simulation
  descrip <- e$sim.realizations %orifnull% 
    if(csr) "simulations of CSR"
    else if(!is.null(simtype)) {
      switch(simtype,
             csr="simulations of CSR",
             rmh=paste("simulations of fitted",
               if(pois) "Poisson" else "Gibbs",
               "model"),
             kppm="simulations of fitted cluster model",
             slrm="simulations of fitted spatial logistic regression model",
             expr="evaluations of user-supplied expression",
             func="evaluations of user-supplied function",
             list="point pattern datasets in user-supplied list",
             funs="columns of user-supplied data",
             "simulated point patterns")
    } else "simulations of fitted model"
  if(!is.null(constraints) && nzchar(constraints))
    descrip <- paste(descrip, constraints)
  #  
  splat("Obtained from", nsim, descrip)
  #
  if(waxlyrical('extras')) {
    if(!is.null(e$dual) && e$dual) 
      splat("Theoretical (i.e. null) mean value of", fname,
            "estimated from a separate set of",
            e$nsim2, "simulations")
    if(!is.null(attr(object, "simfuns")))
      splat("(All", nsim, "simulated function values",
            "are stored in attr(,", dQuote("simfuns"), ") )")
    if(!is.null(attr(object, "simpatterns")))
      splat("(All", nsim, "simulated point patterns",
            "are stored in attr(,", dQuote("simpatterns"), ") )")
  }
  #
  splat("Alternative:", alternative)
  if(V) {
    # nSD envelopes
    splat(switch(alternative,
                 two.sided = "Envelopes",
                 "Critical boundary"),
          "computed as sample mean",
          switch(alternative,
                 two.sided="plus/minus",
                 less="minus",
                 greater="plus"),
          e$nSD, "sample standard deviations")
  } else {
    # critical envelopes
    lo.ord <- if(nr == 1L) "minimum" else paste(ordinal(nr), "smallest")
    hi.ord <- if(nr == 1L) "maximum" else paste(ordinal(nr), "largest")
    if(g) 
      splat(switch(alternative,
                   two.sided = "Envelopes",
                   "Critical boundary"),
            "computed as",
            if(use.theory) "theoretical curve" else "mean of simulations",
            switch(alternative,
                   two.sided="plus/minus",
                   less="minus",
                   greater="plus"),
            hi.ord,
            "simulated value of maximum", 
            switch(alternative,
                   two.sided="absolute",
                   less="negative",
                   greater="positive"),
            "deviation")
    else {
      if(alternative != "less")
        splat("Upper envelope: pointwise", hi.ord, "of simulated curves")
      if(alternative != "greater")
        splat("Lower envelope: pointwise", lo.ord, "of simulated curves")
    }
    symmetric <- (alternative == "two.sided") && !g
    alpha <- if(!symmetric) { nr/(nsim+1) } else { 2 * nr/(nsim+1) }
    splat("Significance level of Monte Carlo test:",
          paste0(if(!symmetric) nr else 2 * nr,
                 "/", nsim+1),
          "=", alpha)
  } 
  splat("Data:", e$Yname)
  return(invisible(NULL))
}
  
envelope.envelope <- function(Y, fun=NULL, ...,
                              transform=NULL, global=FALSE, VARIANCE=FALSE) {

  Yname <- short.deparse(substitute(Y))

  stopifnot(inherits(Y, "envelope"))
  Yorig <- Y

  aargh <- list(...)

  X  <- attr(Y, "datapattern")
  sf <- attr(Y, "simfuns")
  sp <- attr(Y, "simpatterns")
  wt <- attr(Y, "weights")
  einfo <- attr(Y, "einfo")

  csr <- aargh$internal$csr %orifnull% einfo$csr

  if(is.null(fun) && is.null(sf)) {
    # No simulated functions - must compute them from simulated patterns
    if(is.null(sp))
      stop(paste("Cannot compute envelope:",
                 "Y does not contain simulated functions",
                 "(was not generated with savefuns=TRUE)",
                 "and does not contain simulated patterns",
                 "(was not generated with savepatterns=TRUE)"))
    # set default fun=Kest
    fun <- Kest
  }
  
  if(!is.null(fun)) {
    # apply new function 
    # point patterns are required
    if(is.null(sp))
      stop(paste("Object Y does not contain simulated point patterns",
                 "(attribute", dQuote("simpatterns"), ");",
                 "cannot apply a new", sQuote("fun")))
    if(is.null(X))
      stop(paste("Cannot apply a new", sQuote("fun"),
                 "; object Y generated by an older version of spatstat"))
    ## send signal if simulations were CSR
    internal <- aargh$internal
    if(csr) {
        if(is.null(internal)) internal <- list()
        internal$csr <- TRUE
    }
    ## compute new envelope
    result <- do.call(envelope,
                      resolve.defaults(list(Y=quote(X), fun=fun, simulate=sp),
                                       aargh,
                                       list(transform=transform,
                                            global=global,
                                            VARIANCE=VARIANCE,
                                            internal=internal,
                                            Yname=Yname,
                                            nsim=einfo$nsim,
                                            nsim2=einfo$nsim2,
                                            weights=wt),
                                       .StripNull=TRUE))
  } else {
    #' compute new envelope with existing simulated functions
    if(is.null(sf)) 
      stop(paste("Y does not contain a", dQuote("simfuns"), "attribute",
                 "(it was not generated with savefuns=TRUE)"))
    
    if(!is.null(transform)) {
      # Apply transformation to Y and sf
      stopifnot(is.expression(transform))
      ##      cc <- dotexpr.to.call(transform, "Y", "eval.fv")
      cc <- inject.expr("with(Y, .)", transform)
      Y <- eval(cc)
      ##      cc <- dotexpr.to.call(transform, "sf", "eval.fv")
      cc <- inject.expr("with(sf, .)", transform)
      sf <- eval(cc)
    }

    #' catch discrepancy between domains of observed and simulated functions
    if(nrow(sf) != nrow(Y)) {
      rrsim <- sf[[fvnames(sf, ".x")]]
      rrobs <- Y[[fvnames(Y, ".x")]]
      ra <- intersect.ranges(range(rrsim), range(rrobs))
      delta <- min(mean(diff(rrsim)), mean(diff(rrobs)))/2
      oksim <- (rrsim >= ra[1] - delta) & (rrsim <= ra[2] + delta)
      okobs <- (rrobs >= ra[1] - delta) & (rrobs <= ra[2] + delta)
      if(sum(oksim) != sum(okobs))
        stop("Internal error: Unable to reconcile the domains",
             "of the observed and simulated functions", call.=FALSE)
      if(mean(abs(rrsim[oksim] - rrobs[okobs])) > delta)
        stop("Internal error: Unable to reconcile the r values",
             "of the observed and simulated functions", call.=FALSE)
      sf <- sf[oksim, ,drop=FALSE]
      Y  <- Y[okobs,  ,drop=FALSE]
    }
    
    # extract simulated function values
    df <- as.data.frame(sf)
    rname <- fvnames(sf, ".x")
    df <- df[, (names(df) != rname)]

    # interface with 'envelope.matrix'
    etype <- if(global) "global" else if(VARIANCE) "variance" else "pointwise"
    dfm <- as.matrix(df)
    dont.complain.about(dfm)
    result <- do.call(envelope.matrix,
                      resolve.defaults(list(Y=quote(dfm)),
                                       aargh,
                                       list(type=etype,
                                            csr=csr,
                                            funX=Y, 
                                            Yname=Yname,
                                            weights=wt),
                                       .StripNull=TRUE))
  }

  if(!is.null(transform)) {
    # post-process labels
    labl <- attr(result, "labl")
    dnames <- colnames(result)
    dnames <- dnames[dnames %in% fvnames(result, ".")]
    # expand "."
    ud <- as.call(lapply(c("cbind", dnames), as.name))
    dont.complain.about(ud)
    expandtransform <- eval(substitute(substitute(tr, list(.=ud)),
                                       list(tr=transform[[1L]])))
    # compute new labels 
    attr(result, "fname") <- attr(Yorig, "fname")
    mathlabl <- as.character(fvlegend(result, expandtransform))
    # match labels to columns
    evars <- all.vars(expandtransform)
    used.dotnames <- evars[evars %in% dnames]
    mathmap <- match(colnames(result), used.dotnames)
    okmath <- !is.na(mathmap)
    # update appropriate labels
    labl[okmath] <- mathlabl[mathmap[okmath]]
    attr(result, "labl") <- labl
  }
  
  # Tack on envelope info
  copyacross <- c("Yname", "csr.theo", "use.theory", "simtype", "constraints")
  attr(result, "einfo")[copyacross] <- attr(Yorig, "einfo")[copyacross]
  attr(result, "einfo")$csr <- csr
  # Save data
  
  return(result)
}

pool.envelope <- local({

  pool.envelope <- function(..., savefuns=FALSE, savepatterns=FALSE) {
    Yname <- short.deparse(sys.call())
    if(nchar(Yname) > 60) Yname <- paste(substr(Yname, 1L, 40L), "[..]")
    Elist <- unname(list(...))
    nE <-  length(Elist)
    if(nE == 0) return(NULL)
    #' ........ validate envelopes .....................
    #' All arguments must be envelopes
    notenv <- !unlist(lapply(Elist, inherits, what="envelope"))
    if(any(notenv)) {
      n <- sum(notenv)
      why <- paste(ngettext(n, "Argument", "Arguments"),
                   commasep(which(notenv)),
                   ngettext(n, "does not", "do not"),
                   "belong to the class",
                   dQuote("envelope"))
      stop(why)
    }
    E1 <- Elist[[1L]]
    ## Only one envelope?
    if(nE == 1)
      return(E1)
    ## envelopes must be compatible
    ok <- do.call(compatible, Elist)
    if(!ok)
      stop("Envelopes are not compatible")
    ## find name of function argument
    argname <- fvnames(E1, ".x")
    arg.desc <- attr(E1, "desc")[match(argname, colnames(E1))]
    ## ... reconcile parameters in different envelopes .......
    eilist <- lapply(Elist, attr, which="einfo")
    nrank  <- resolveEinfo(eilist, "nrank", 1)
    global    <- resolveEinfo(eilist, "global",   FALSE)
    ginterval    <- resolveEinfo(eilist, "ginterval", NULL, atomic=FALSE)
    VARIANCE  <- resolveEinfo(eilist, "VARIANCE", FALSE)
    alternative      <- resolveEinfo(eilist, "alternative", FALSE)
    scale <- resolveEinfo(eilist, "scale", NULL, atomic=FALSE)
    clamp <- resolveEinfo(eilist, "clamp", FALSE)
    resolveEinfo(eilist, "simtype",  "funs",
                 "Envelopes were generated using different types of simulation")
    resolveEinfo(eilist, "constraints",  "",
            "Envelopes were generated using different types of conditioning")
    resolveEinfo(eilist, "csr.theo", FALSE, NULL)
    csr         <- resolveEinfo(eilist, "csr", FALSE, NULL)
    use.weights <- resolveEinfo(eilist, "use.weights" , FALSE,
     "Weights were used in some, but not all, envelopes: they will be ignored")
    use.theory <- resolveEinfo(eilist, "use.theory", csr, NULL)
    ##
    weights <-
      if(use.weights) unlist(lapply(Elist, attr, which="weights")) else NULL
    type <- if(global) "global" else if(VARIANCE) "variance" else "pointwise"
    
    ## ........ validate saved functions .....................
    if(savefuns || !VARIANCE) {
      ## Individual simulated functions are required
      SFlist <- lapply(Elist, attr, which="simfuns")
      isnul <- unlist(lapply(SFlist, is.null))
      if(any(isnul)) {
        n <- sum(isnul)
        comply <- if(!VARIANCE) "compute the envelope:" else
                  "save the simulated functions:"
        why <- paste("Cannot", comply,
                     ngettext(n, "argument", "arguments"),
                     commasep(which(isnul)),
                     ngettext(n, "does not", "do not"),
                     "contain a", dQuote("simfuns"), "attribute",
                     "(not generated with savefuns=TRUE)")
        stop(why)
      }
      ## Simulated functions must be the same function
      fnames <- unique(lapply(SFlist, attr, which="fname"))
      if(length(fnames) > 1L) {
        fnames <- unlist(lapply(fnames, flatfname))
        stop(paste("Envelope objects contain values",
                   "of different functions:",
                   commasep(sQuote(fnames))))
      }
      ## vectors of r values must be identical
      rlist <- lapply(SFlist, getargvals)
      argvals <- rlist[[1L]]
      samer <- unlist(lapply(rlist, identical, y=argvals))
      if(!all(samer))
        stop(paste("Simulated function values are not compatible",
                   "(different values of function argument)"))
      ## Extract function values and assemble into one matrix
      matlist <- lapply(SFlist, getdotvals)
      SFmatrix <- do.call(cbind, matlist)
    }
    ## compute pooled envelope
    switch(type,
           pointwise = {
             result <- envelope(SFmatrix, funX=E1,
                                type=type, alternative=alternative,
                                clamp=clamp,
                                nrank=nrank,
                                csr=csr, use.theory=use.theory,
                                Yname=Yname, weights=weights,
                                savefuns=savefuns)
           },
           global = {
             simfunmatrix <- if(is.null(ginterval)) SFmatrix else {
               ## savefuns have not yet been clipped to ginterval
               ## while envelope data have been clipped.
               domain <- (argvals >= ginterval[1L]) & (argvals <= ginterval[2L])
               SFmatrix[domain, , drop=FALSE]
             }
             result <- envelope(simfunmatrix, funX=E1,
                                type=type, alternative=alternative,
                                scale=scale, clamp=clamp,
                                csr=csr, use.theory=use.theory,
                                nrank=nrank,
                                ginterval=ginterval,
                                Yname=Yname, weights=weights,
                                savefuns=savefuns)
           },
           variance = {
             ## Pool sample means and variances
             nsims <- unlist(lapply(eilist, getElement, name="nsim"))
             mmeans <- lapply(Elist, getElement, name="mmean")
             vars   <- lapply(Elist, getElement, name="var")
             mmeans <- matrix(unlist(mmeans), ncol=nE)
             vars   <- matrix(unlist(vars),   ncol=nE)
             if(!use.weights) {
               w.mean <- nsims
               d.mean <- sum(nsims)
               w.var  <- nsims - 1
               d.var  <- sum(nsims) - 1
             } else {
               weightlist <- lapply(Elist, attr, which="weights")
               w.mean <- unlist(lapply(weightlist, sum))
               d.mean <- sum(w.mean)
               ssw <- unlist(lapply(weightlist, meansqfrac))
               ##  meansqfrac :  function(x) {sum((x/sum(x))^2)}))
               w.var  <- w.mean * (1 - ssw)
               d.var <-  d.mean * (1 - sum(ssw))
             }
             poolmmean <- as.numeric(mmeans %*% matrix(w.mean/d.mean, ncol=1L))
             within <- vars %*% matrix(w.var, ncol=1L)
             between <- ((mmeans - poolmmean[])^2) %*% matrix(w.mean, ncol=1L)
             poolvar <- as.numeric((within + between)/d.var)
             ## feed precomputed data to envelope.matrix
             pc <- list(Ef=poolmmean[],
                        varf=poolvar[])
             nsim <- sum(nsims)
             result <- envelope.matrix(NULL, funX=E1,
                                       type=type, alternative=alternative,
                                       csr=csr, Yname=Yname,
                                       weights=weights,
                                       savefuns=savefuns,
                                       nsim=nsim,
                                       precomputed=pc)
           })
  
    ## Copy envelope info that is not handled by envelope.matrix
    copyacross <- c("Yname", "csr.theo", "use.theory", "simtype", "constraints")
    attr(result, "einfo")[copyacross] <- attr(E1, "einfo")[copyacross]
  
    ## ..............saved patterns .....................
    if(savepatterns) {
      SPlist <- lapply(Elist, attr, which="simpatterns")
      isnul <- unlist(lapply(SPlist, is.null))
      if(any(isnul)) {
        n <- sum(isnul)
        why <- paste("Cannot save the simulated patterns:",
                     ngettext(n, "argument", "arguments"),
                     commasep(which(isnul)),
                     ngettext(n, "does not", "do not"),
                     "contain a", dQuote("simpatterns"), "attribute",
                     "(not generated with savepatterns=TRUE)")
        warning(why)
      } else {
        attr(result, "simpatterns") <- Reduce(append, SPlist)
      }
    }
    ## ..............saved summary functions ................
    if(savefuns) {
      alldata <- cbind(argvals, SFmatrix)
      Nsim <- ncol(SFmatrix)
      simnames <- paste0("sim", 1:Nsim)
      colnames(alldata) <- c(argname, simnames)
      alldata <- as.data.frame(alldata)
      SFtemplate <- SFlist[[1L]]
      SimFuns <- fv(alldata,
                    argu=argname,
                    ylab=attr(SFtemplate, "ylab"),
                    valu="sim1",
                    fmla=paste(". ~", argname),
                    alim=attr(SFtemplate, "alim"),
                    labl=names(alldata),
                    desc=c(arg.desc,
                      paste("Simulation ", 1:Nsim, sep="")),
                    fname=attr(SFtemplate, "fname"),
                    yexp=attr(SFtemplate, "yexp"),
                    unitname=unitname(SFtemplate))
      fvnames(SimFuns, ".") <- simnames
      attr(result, "simfuns") <- SimFuns
    } 
  
    dotnames   <- lapply(Elist, fvnames, a=".")
    dn <- dotnames[[1L]]
    if(all(unlist(lapply(dotnames, identical, y=dn))))
      fvnames(result, ".") <- dn
    
    shadenames <- lapply(Elist, fvnames, a=".s")
    sh <- shadenames[[1L]]
    if(all(unlist(lapply(shadenames, identical, y=sh))))
      fvnames(result, ".s") <- sh
  
    return(result)
  }

  getargvals <- function(z) { as.matrix(z)[, fvnames(z, ".x")] }
  
  getdotvals <- function(z) { as.matrix(z)[, fvnames(z, "."), drop=FALSE] }

  meansqfrac <- function(x) {sum((x/sum(x))^2)} 

  pool.envelope
})

# resolve matching entries in different envelope objects
#   x is a list of envelope info objects

resolveEinfo <- function(x, what, fallback, warn, atomic=TRUE) {
  if(atomic) {
    y <- unique(unlist(lapply(x, getElement, name=what)))
    if(length(y) == 1L)
      return(y)
  } else {
    y <- unique(lapply(x, getElement, name=what))
    if(length(y) == 1L)
      return(y[[1L]])
  }
  if(missing(warn))
    warn <- paste("Envelopes were generated using different values",
                  "of argument", paste(sQuote(what), ";", sep=""),
                  "reverting to default value")
  if(!is.null(warn))
    warning(warn, call.=FALSE)
  return(fallback)
}

