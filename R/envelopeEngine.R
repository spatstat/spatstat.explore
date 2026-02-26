#'
#'    envelopeEngine.R
#'
#'    Main internal code for simulating and computing envelopes
#'
#'    envelopeEngine()    performs simulations and calculates summary function
#'    envelope.matrix()   calculates envelope from function values
#'
#'    $Revision: 1.7 $ $Date: 2026/02/26 09:32:42 $
#' 
#' ............. envelopeEngine() ..................................
#'
#' X is the data point pattern, which could be ppp, pp3, ppx etc
#' X determines the class of pattern expected from the simulations

envelopeEngine <-
  function(X, fun, simul,
           nsim=99, nrank=1, ..., funargs=list(), funYargs=funargs,
           verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           theoryfun=NULL, theory.adjective=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE,
           saveresultof=NULL,
           weights=NULL,
           nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL,
           maxnerr=nsim,
           rejectNA=FALSE,
           silent=FALSE,
           maxerr.action=c("fatal", "warn", "null"),
           internal=NULL, cl=NULL,
           envir.user=envir.user,
           expected.arg="r",
           do.pwrong=FALSE,
           foreignclass=NULL,
           collectrubbish=FALSE)
{
  #
  envir.here <- sys.frame(sys.nframe())

  alternative <- match.arg(alternative)
  maxerr.action <- match.arg(maxerr.action)  

  foreignclass <- as.character(foreignclass)
  if(length(foreignclass) != 0 && clipdata) {
    warning(paste("Ignoring clipdata=TRUE:",
                  "I don't know how to clip objects of class",
                  sQuote(paste(foreignclass, collapse=","))))
    clipdata <- FALSE
  }
  
  # ----------------------------------------------------------
  # Determine Simulation
  # ----------------------------------------------------------

  check.1.integer(nsim)
  stopifnot(nsim > 1)

  # Identify class of patterns to be simulated, from data pattern X
  Xclass <- if(is.ppp(X)) "ppp" else
            if(is.pp3(X)) "pp3" else
            if(is.ppx(X)) "ppx" else
            if(inherits(X, foreignclass)) foreignclass else
            stop("Unrecognised class of point pattern")
  Xobjectname <- paste("point pattern of class", sQuote(Xclass))

  # Option to use weighted average
  if(use.weights <- !is.null(weights)) {
    # weight can be either a numeric vector or a function
    if(is.numeric(weights)) {
      compute.weights <- FALSE
      weightfun <- NULL
    } else if(is.function(weights)) {
      compute.weights <- TRUE
      weightfun <- weights
      weights <- NULL  
    } else stop("weights should be either a function or a numeric vector")
  } else compute.weights <- FALSE
    
  # Undocumented option to generate patterns only.
  patterns.only <- identical(internal$eject, "patterns")

  # Undocumented option to evaluate 'something' for each simulation
  if(savevalues <- !is.null(saveresultof)) {
    stopifnot(is.function(saveresultof))
    SavedValues <- list()
    ## might be a function of the pattern only, or both pattern and summary fun
    result.depends.both <- (length(formals(saveresultof)) >= 2)
  }

  # Identify type of simulation from argument 'simul'
  if(inherits(simul, "simulrecipe")) {
    # ..................................................
    # simulation recipe is given
    simtype <- simul$type
    simexpr <- simul$expr
    envir   <- simul$envir
    csr     <- simul$csr
    pois    <- simul$pois
    constraints <- simul$constraints
  } else {
    # ...................................................
    # simulation is specified by argument `simulate' to envelope()
    simulate <- simul
    # which should be an expression, or a list of point patterns,
    # or an envelope object, or a function to be applied to the data
    csr <- FALSE
    # override
    if(!is.null(icsr <- internal$csr)) csr <- icsr
    pois <- csr
    constraints <- ""
#    model <- NULL
    if(inherits(simulate, "envelope")) {
      # envelope object: see if it contains stored point patterns
      simpat <- attr(simulate, "simpatterns")
      if(!is.null(simpat))
        simulate <- simpat
      else
        stop(paste("The argument", sQuote("simulate"),
                   "is an envelope object but does not contain",
                   "any saved point patterns."))
    }
    if(is.expression(simulate)) {
      ## The user-supplied expression 'simulate' will be evaluated repeatedly
      simtype <- "expr"
      simexpr <- simulate
      envir <- envir.user
    } else if(is.function(simulate)) {
      ## User-supplied function 'simulate' will be repeatedly evaluated on X
      simtype <- "func"
      simexpr <- expression(simulate(X))
      envir <- envir.here
    } else if(is.list(simulate) &&
              all(sapply(simulate, inherits, what=Xclass))) {
      #' The user-supplied list of point patterns will be used
      simtype <- "list"
      SimDataList <- simulate
      #' expression that will be evaluated
      simexpr <- expression(SimDataList[[i+nerr]])
      dont.complain.about(SimDataList)
      envir <- envir.here
      #' ensure that `i' is defined
      i <- 1L
      nerr <- 0L
      maxnerr <- min(length(SimDataList)-nsim, maxnerr)
      #' any messages?
      if(!is.null(mess <- attr(simulate, "internal"))) {
        # determine whether these point patterns are realisations of CSR
        csr <- isTRUE(mess$csr)
      }
    } else if(is.list(simulate) &&
              all(sapply(simulate, is.list)) &&
              all(lengths(simulate) == 1) &&
              all(sapply((elements <- lapply(simulate, "[[", i=1)), inherits, what=Xclass))) {
      #' malformed argument: list(list(ppp), list(ppp), ....) 
      SimDataList <- elements
      simtype <- "list"
      #' expression that will be evaluated
      simexpr <- expression(SimDataList[[i+nerr]])
      dont.complain.about(SimDataList)
      envir <- envir.here
      #' ensure that `i' is defined
      i <- 1L
      nerr <- 0L
      maxnerr <- min(length(SimDataList)-nsim, maxnerr)
      #' any messages?
      if(!is.null(mess <- attr(simulate, "internal"))) {
        # determine whether these point patterns are realisations of CSR
        csr <- isTRUE(mess$csr) 
      }
    } else stop(paste(sQuote("simulate"),
                      "should be an expression,",
                      "or a list of point patterns of the same kind as X"))
  }
  # -------------------------------------------------------------------
  # Determine clipping window
  # ------------------------------------------------------------------

  if(clipdata) {
    # Generate one realisation
    Xsim <- eval(simexpr, envir=envir)
    if(!inherits(Xsim, Xclass))
      switch(simtype,
             csr=stop(paste("Internal error:", Xobjectname, "not generated")),
             rmh=stop(paste("Internal error: rmh did not return an",
               Xobjectname)),
             kppm=stop(paste("Internal error: simulate.kppm did not return an",
               Xobjectname)),
             slrm=stop(paste("Internal error: simulate.slrm did not return an",
               Xobjectname)),
             expr=stop(paste("Evaluating the expression", sQuote("simulate"),
               "did not yield an", Xobjectname)),
             func=stop(paste("Evaluating the function", sQuote("simulate"),
               "did not yield an", Xobjectname)),
             list=stop(paste("Internal error: list entry was not an",
               Xobjectname)),
             stop(paste("Internal error:", Xobjectname, "not generated"))
             )
    # Extract window
    clipwin <- Xsim$window
    if(!is.subset.owin(clipwin, X$window))
      warning("Window containing simulated patterns is not a subset of data window")
  }
  
  # ------------------------------------------------------------------
  # Summary function to be applied 
  # ------------------------------------------------------------------

  if(is.null(fun))
    stop("Internal error: fun is NULL")

  # Name of function, for error messages
  fname <- if(is.name(substitute(fun))) short.deparse(substitute(fun)) else
  if(is.character(fun)) fun else "fun"
  fname <- sQuote(fname)

  # R function to apply
  if(is.character(fun)) {
    gotfun <- try(get(fun, mode="function"))
    if(inherits(gotfun, "try-error"))
      stop(paste("Could not find a function named", sQuote(fun)))
    fun <- gotfun
  } else if(!is.function(fun)) 
    stop(paste("unrecognised format for function", fname))
  fargs <- names(formals(fun))
  if(!any(c(expected.arg, "...") %in% fargs))
    stop(paste(fname, "should have",
               ngettext(length(expected.arg), "an argument", "arguments"),
               "named", commasep(sQuote(expected.arg)),
               "or a", sQuote("..."), "argument"))
  usecorrection <- any(c("correction", "...") %in% fargs)
  usezerocor <- any(c("zerocor", "...") %in% fargs)
  
  # ---------------------------------------------------------------------
  # validate other arguments
  if((nrank %% 1) != 0)
    stop("nrank must be an integer")
  if((nsim %% 1) != 0)
    stop("nsim must be an integer")
  stopifnot(nrank > 0 && nrank < nsim/2)

  arg.given <- any(expected.arg %in% names(list(...)))

  if(tran <- !is.null(transform)) {
    stopifnot(is.expression(transform))
    # prepare expressions to be evaluated each time 
    transform.funX    <- inject.expr("with(funX,.)",    transform)
    transform.funXsim <- inject.expr("with(funXsim,.)", transform)
    # .... old code using 'eval.fv' ......
    # transform.funX <- dotexpr.to.call(transform, "funX", "eval.fv")
    # transform.funXsim <- dotexpr.to.call(transform, "funXsim", "eval.fv")
    # 'transform.funX' and 'transform.funXsim' are unevaluated calls to eval.fv
  }
  if(!is.null(ginterval)) 
    stopifnot(is.numeric(ginterval) && length(ginterval) == 2)
    
  # ---------------------------------------------------------------------
  # Evaluate function for data pattern X
  # ------------------------------------------------------------------
  Xarg <- if(!clipdata) X else X[clipwin]
  corrx <- if(usecorrection) list(correction="best") else NULL
  corrz <- if(usezerocor) list(zerocor="best") else NULL
  dont.complain.about(Xarg)
  funX <- do.call(fun,
                  resolve.defaults(list(quote(Xarg)),
                                   list(...),
                                   funYargs,
                                   corrx,
                                   corrz))
                                     
  if(!inherits(funX, "fv"))
    stop(paste("The function", fname,
               "must return an object of class", sQuote("fv")))

  ## catch 'conservation' parameters
  conserveargs <- attr(funX, "conserve")
  if(!is.null(conserveargs) && !any(c("conserve", "...") %in% fargs))
    stop(paste("In this usage, the function", fname,
               "should have an argument named 'conserve' or '...'"))

  ## warn about 'dangerous' arguments
  if(!is.null(dang <- attr(funX, "dangerous")) &&
     any(uhoh <- dang %in% names(list(...)))) {
    nuh <- sum(uhoh)
    warning(paste("Envelope may be invalid;",
                  ngettext(nuh, "argument", "arguments"),
                  commasep(sQuote(dang[uhoh])),
                  ngettext(nuh, "appears", "appear"),
                  "to have been fixed."),
            call.=FALSE)
  }
  
  argname <- fvnames(funX, ".x")
  valname <- fvnames(funX, ".y")
  has.theo <- "theo" %in% fvnames(funX, "*")

  if(compute.theo <- is.function(theoryfun)) {
    ## Theoretical value is provided by a function
    theorydesc <- pasteN(theory.adjective, "theoretical value of %s")
    if(has.theo) {
      ## overwrite column values
      rvals <- funX[[argname]]
      funX[["theo"]] <- theoryfun(rvals)
      funX <- tweak.fv.entry(funX, "theo", new.desc=theorydesc)
    } else {
      ## create new column
      funX <- bind.fv(funX, theoryfun, labl="theo", desc=theorydesc)
      has.theo <- TRUE
    }
  }
  csr.theo <- csr && has.theo
  use.theory <- compute.theo || (has.theo && (use.theory %orifnull% csr))
  
  if(tran) {
    # extract only the recommended value
    if(use.theory) 
      funX <- funX[, c(argname, valname, "theo")]
    else
      funX <- funX[, c(argname, valname)]
    ## save name of function before transformation
    fname.orig <- attr(funX, "fname")
    ## apply the transformation to it
    funX <- eval(transform.funX)
  } else {
    fname.orig <- NULL
  }
    
  argvals <- funX[[argname]]
#  fX    <- funX[[valname]]

  arg.desc <- attr(funX, "desc")[match(argname, colnames(funX))]
  
  # default domain over which to maximise
  alim <- attr(funX, "alim")
  if(global && is.null(ginterval))
    ginterval <- if(arg.given || is.null(alim)) range(argvals) else alim
  
  #--------------------------------------------------------------------
  # Determine number of simulations
  # ------------------------------------------------------------------
  #
  ## determine whether dual simulations are required
  ## (one set of simulations to calculate the theoretical mean,
  ##  another independent set of simulations to obtain the critical point.)
  dual <- (global && !use.theory && !VARIANCE)
  if(dual) {
    check.1.integer(nsim2)
    stopifnot(nsim2 >= 1)
  }
  Nsim <- if(!dual) nsim else (nsim + nsim2)

  # if taking data from a list of point patterns,
  # check there are enough of them
  if(simtype == "list" && Nsim > length(SimDataList))
    stop(paste("Number of simulations",
               paren(if(!dual)
                     paste(nsim) else
                     paste(nsim, "+", nsim2, "=", Nsim)
                     ),
               "exceeds number of point pattern datasets supplied"))

  # Undocumented secret exit
  # ------------------------------------------------------------------
  if(patterns.only) {
    # generate simulated realisations and return only these patterns
    if(verbose) {
      action <- if(simtype == "list") "Extracting" else "Generating"
      descrip <- switch(simtype,
                        csr = "simulations of CSR",
                        rmh = paste("simulated realisations of fitted",
                          if(pois) "Poisson" else "Gibbs",
                          "model"),
                        kppm = "simulated realisations of fitted cluster model",
                        slrm = "simulated realisations of spatial logistic regression model",
                        expr = "simulations by evaluating expression",
                        func = "simulations by evaluating function",
                        list = "point patterns from list",
                        "simulated realisations")
      if(!is.null(constraints) && nzchar(constraints))
        descrip <- paste(descrip, constraints)
      explan <- if(dual) paren(paste(nsim2, "to estimate the mean and",
                                     nsim, "to calculate envelopes")) else ""
      splat(action, Nsim, descrip, explan, "...")
    }
    XsimList <- list()
  # start simulation loop
    sstate <- list() 
    for(i in 1:Nsim) {
      if(verbose) sstate <- progressreport(i, Nsim, state=sstate)
      Xsim <- eval(simexpr, envir=envir)
      if(!inherits(Xsim, Xclass))
        switch(simtype,
               csr={
                 stop(paste("Internal error:", Xobjectname, "not generated"))
               },
               rmh={
                 stop(paste("Internal error: rmh did not return an",
                            Xobjectname))
               },
               kppm={
                 stop(paste("Internal error: simulate.kppm did not return an",
                            Xobjectname))
               },
               slrm={
                 stop(paste("Internal error: simulate.slrm did not return an",
                            Xobjectname))
               },
               expr={
                 stop(paste("Evaluating the expression", sQuote("simulate"),
                            "did not yield an", Xobjectname))
               },
               func={
                 stop(paste("Evaluating the function", sQuote("simulate"),
                            "did not yield an", Xobjectname))
               },
               list={
                 stop(paste("Internal error: list entry was not an",
                            Xobjectname))
               },
               stop(paste("Internal error:", Xobjectname, "not generated"))
               )
      XsimList[[i]] <- Xsim
    }
    if(verbose) {
      cat(paste("Done.\n"))
      flush.console()
    }
    attr(XsimList, "internal") <- list(csr=csr)
    return(XsimList)
  }
  
  # capture main decision parameters
  envelopeInfo <-  list(call=cl,
                        Yname=Yname,
                        valname=valname,
                        csr=csr,
                        csr.theo=csr.theo,
                        use.theory=use.theory,
                        theoryfun=theoryfun,
                        theory.adjective=theory.adjective,
                        pois=pois,
                        simtype=simtype,
                        constraints=constraints,
                        nrank=nrank,
                        nsim=nsim,
                        Nsim=Nsim,
                        global=global,
                        ginterval=ginterval,
                        dual=dual,
                        nsim2=nsim2,
                        VARIANCE=VARIANCE,
                        nSD=nSD,
                        alternative=alternative,
                        scale=scale,
                        clamp=clamp,
                        use.weights=use.weights,
                        do.pwrong=do.pwrong)

  # ----------------------------------------
  ######### SIMULATE #######################
  # ----------------------------------------

  if(verbose) {
    action <- if(simtype == "list") "Extracting" else "Generating"
    descrip <- switch(simtype,
                      csr = "simulations of CSR",
                      rmh = paste("simulated realisations of fitted",
                        if(pois) "Poisson" else "Gibbs",
                        "model"),
                      kppm = "simulated realisations of fitted cluster model",
                      slrm = "simulated realisations of fitted spatial logistic regression model",
                      expr = "simulations by evaluating expression",
                      func = "simulations by evaluating function",
                      list = "point patterns from list",
                      "simulated patterns")
    if(!is.null(constraints) && nzchar(constraints))
      descrip <- paste(descrip, constraints)
    explan <- if(dual) paren(paste(nsim2, "to estimate the mean and",
                                   nsim, "to calculate envelopes")) else ""
    splat(action, Nsim, descrip, explan, "...")
  }
  # determine whether simulated point patterns should be saved
  catchpatterns <- savepatterns && simtype != "list"
  Caughtpatterns <- list()
  # allocate space for computed function values
  nargvals <- length(argvals)
  simvals <- matrix(, nrow=nargvals, ncol=Nsim)
  # allocate space for weights to be computed
  if(compute.weights)
    weights <- numeric(Nsim)
  
  # inferred values of function argument 'r' or equivalent parameters
  if(identical(expected.arg, "r") && identical(argname, "r")) {
    # Kest, etc
    inferred.r.args <- list(r=argvals)
  } else if(identical(expected.arg, c("rmax", "nrval"))) {
    # K3est, etc
    inferred.r.args <- list(rmax=max(argvals), nrval=length(argvals))
  } else if(any(c(argname, "...") %in% fargs)) {
    ## assume it accepts the vector of argument values
    inferred.r.args <- structure(list(argvals), names=argname)
  } else {
    stop(paste("Don't know how to infer values of",
               if(length(expected.arg))
                 commasep(sQuote(expected.arg)) else "function argument"))
  }
    
  # arguments for function when applied to simulated patterns
  funargs <-
    resolve.defaults(funargs,
                     inferred.r.args,
                     list(...),
                     conserveargs,
                     if(usecorrection) list(correction="best") else NULL,
                     if(usezerocor) list(zerocor="best") else NULL)

  # reject simulated pattern if function values are all NA (etc)
  rejectNA <- isTRUE(rejectNA)
  
  # start simulation loop
  nerr <- 0
  gaveup <- FALSE
  if(verbose) pstate <- list()
  for(i in 1:Nsim) {
    ## safely generate a random pattern and apply function
    success <- FALSE
    while(!success && !gaveup) {
      Xsim <- eval(simexpr, envir=envir)
      ## check valid point pattern
      if(!inherits(Xsim, Xclass))
        switch(simtype,
               csr=stop(paste("Internal error:", Xobjectname, "not generated")),
               rmh=stop(paste("Internal error: rmh did not return an",
                 Xobjectname)),
               kppm=stop(paste("Internal error:",
                 "simulate.kppm did not return an",
                 Xobjectname)),
               slrm=stop(paste("Internal error:",
                 "simulate.slrm did not return an",
                 Xobjectname)),
               expr=stop(paste("Evaluating the expression", sQuote("simulate"),
                 "did not yield an", Xobjectname)),
               func=stop(paste("Evaluating the function", sQuote("simulate"),
                 "did not yield an", Xobjectname)),
               list=stop(paste("Internal error: list entry was not an",
                 Xobjectname)),
               stop(paste("Internal error:", Xobjectname, "not generated"))
               )
      if(catchpatterns)
        Caughtpatterns[[i]] <- Xsim
      if(savevalues && !result.depends.both)
        SavedValues[[i]] <- saveresultof(Xsim)
      if(compute.weights) {
        wti <- weightfun(Xsim)
        if(!is.numeric(wti))
          stop("weightfun did not return a numeric value")
        if(length(wti) != 1L)
          stop("weightfun should return a single numeric value")
        weights[i] <- wti
      }
      ## apply function safely
      funXsim <- try(do.call(fun, c(list(Xsim), funargs)), silent=silent)

      success <-
        !inherits(funXsim, "try-error") &&
        inherits(funXsim, "fv") &&
        (!rejectNA || any(is.finite(funXsim[[valname]])))

      if(!success) {
        #' error in computing summary function
        nerr <- nerr + 1L 
        if(nerr > maxnerr) {
          gaveup <- TRUE
          errtype <- if(rejectNA) "fatal errors or NA function values"
          if(simtype == "list") {
            whinge <- paste("Exceeded maximum possible number of errors",
                          "when evaluating summary function:",
                          length(SimDataList), "patterns provided,",
                          nsim, "patterns required,",
                          nerr, ngettext(nerr, "pattern", "pattern"),
                          "rejected due to", errtype)
          } else {
            whinge <- paste("Exceeded maximum permissible number of",
                            errtype,
                            paren(paste("maxnerr =", maxnerr)),
                            "when evaluating summary function",
                            "for simulated point patterns")
          }
          switch(maxerr.action,
                 fatal = stop(whinge, call.=FALSE),
                 warn  = warning(whinge, call.=FALSE),
                 null  = {})
        } else if(!silent) cat("[retrying]\n")
      }

      #' ..... end of while(!success) ................
    }
    
    if(gaveup) break; # exit loop now
    
    ## sanity checks
    if(i == 1L) {
      if(!inherits(funXsim, "fv"))
        stop(paste("When applied to a simulated pattern, the function",
                   fname, "did not return an object of class",
                   sQuote("fv")))
      argname.sim <- fvnames(funXsim, ".x")
      valname.sim <- fvnames(funXsim, ".y")
      if(argname.sim != argname)
        stop(paste("The objects returned by", fname,
                   "when applied to a simulated pattern",
                   "and to the data pattern",
                   "are incompatible. They have different argument names",
                   sQuote(argname.sim), "and", sQuote(argname), 
                   "respectively"))
      if(valname.sim != valname)
        stop(paste("When", fname, "is applied to a simulated pattern",
                   "it provides an estimate named", sQuote(valname.sim), 
                   "whereas the estimate for the data pattern is named",
                   sQuote(valname),
                   ". Try using the argument", sQuote("correction"),
                   "to make them compatible"))
      rfunX    <- with(funX,    .x)
      rfunXsim <- with(funXsim, .x)
      if(!identical(rfunX, rfunXsim))
        stop(paste("When", fname, "is applied to a simulated pattern,",
                   "the values of the argument", sQuote(argname.sim),
                   "are different from those used for the data."))
    }

    if(compute.theo) {
      ## Theoretical value is provided by a function
      if("theo" %in% fvnames(funXsim, "*")) {
        ## overwrite existing column values
        funXsim[["theo"]] <- theoryfun(rfunXsim)
        funXsim <- tweak.fv.entry(funXsim, "theo", new.desc=theorydesc)
      } else {
        ## create new column
        funXsim <- bind.fv(funXsim, theoryfun,
                           labl="theo", desc=theorydesc)
      }
    }
    
    if(savevalues && result.depends.both)
      SavedValues[[i]] <- saveresultof(Xsim, funXsim)

    if(tran) {
      # extract only the recommended value
      if(use.theory) 
        funXsim <- funXsim[, c(argname, valname, "theo")]
      else
        funXsim <- funXsim[, c(argname, valname)]
      # apply the transformation to it
      funXsim <- eval(transform.funXsim)
    }

    # extract the values for simulation i
    simvals.i <- funXsim[[valname]]
    if(length(simvals.i) != nargvals)
      stop("Vectors of function values have incompatible lengths")
      
    simvals[ , i] <- funXsim[[valname]]
    if(verbose)
      pstate <- progressreport(i, Nsim, state=pstate)
    
    if(collectrubbish) {
      rm(Xsim)
      rm(funXsim)
      gc()
    }
  }
  ##  end simulation loop
  
  if(verbose) {
    cat("\nDone.\n")
    flush.console()
  }

  # ...........................................................
  # save functions and/or patterns if so commanded

  if(!gaveup) {
    if(savefuns) {
      alldata <- cbind(argvals, simvals)
      simnames <- paste("sim", 1:Nsim, sep="")
      colnames(alldata) <- c(argname, simnames)
      alldata <- as.data.frame(alldata)
      SimFuns <- fv(alldata,
                    argu=argname,
                    ylab=attr(funX, "ylab"),
                    valu="sim1",
                    fmla= paste(". ~", argname),
                    alim=attr(funX, "alim"),
                    labl=names(alldata),
                    desc=c(arg.desc,
                           paste("Simulation ", 1:Nsim, sep="")),
                    fname=attr(funX, "fname"),
                    yexp=attr(funX, "yexp"),
                    unitname=unitname(funX))
      fvnames(SimFuns, ".") <- simnames
    } 
    if(savepatterns)
      SimPats <- if(simtype == "list") SimDataList else Caughtpatterns
  }
  
  ######### COMPUTE ENVELOPES #######################

  etype <- if(global) "global" else if(VARIANCE) "variance" else "pointwise"
  if(dual) {
    jsim <- 1:nsim
    jsim.mean <- nsim + 1:nsim2
  } else {
    jsim <- jsim.mean <- NULL
  }

  result <- envelope.matrix(simvals, funX=funX,
                            jsim=jsim, jsim.mean=jsim.mean,
                            type=etype, alternative=alternative,
                            scale=scale, clamp=clamp,
                            csr=csr, use.theory=use.theory,
                            theory.adjective=theory.adjective,
                            nrank=nrank, ginterval=ginterval, nSD=nSD,
                            fname.orig=fname.orig, transform=transform,
                            Yname=Yname, do.pwrong=do.pwrong,
                            weights=weights, gaveup=gaveup)

  ## tack on envelope information
  attr(result, "einfo") <- resolve.defaults(envelopeInfo,
                                            attr(result, "einfo"))

  if(!gaveup) {
    ## tack on functions and/or patterns if so commanded   
    if(savefuns)
      attr(result, "simfuns") <- SimFuns
    if(savepatterns) {
      attr(result, "simpatterns") <- SimPats
      attr(result, "datapattern") <- X
    }
    ## undocumented - tack on values of some other quantity
    if(savevalues) {
      attr(result, "simvalues") <- SavedValues
      attr(result, "datavalue") <-
        if(result.depends.both) saveresultof(X, funX) else saveresultof(X)
    }
  }

  ## save function weights 
  if(use.weights)
    attr(result, "weights") <- weights
  return(result)
}


#'
#' ................ envelope.matrix ..........................
#'
#'     core functionality to compute envelope values from function values
#'
#'     theory   = funX[["theo"]]
#'     observed = fX

envelope.matrix <- function(Y, ...,
                            argvals=rvals, rvals=NULL, ## rvals is old name
                            observed=NULL, theory=NULL,
                            funX=NULL,
                            nsim=NULL, nsim2=NULL,
                            jsim=NULL, jsim.mean=NULL,
                            type=c("pointwise", "global", "variance"),
                            alternative=c("two.sided", "less", "greater"),
                            scale = NULL, clamp=FALSE,
                            csr=FALSE, use.theory = csr, 
                            nrank=1, ginterval=NULL, nSD=2,
                            savefuns=FALSE,
                            check=TRUE,
                            Yname=NULL,
                            argname=NULL,
                            arg.desc=NULL,
                            theory.adjective=NULL,
                            fname.orig=NULL,
                            transform=NULL,
                            do.pwrong=FALSE,
                            weights=NULL,
                            precomputed=NULL,
                            gaveup=FALSE) {
  if(is.null(Yname))
    Yname <- short.deparse(substitute(Y))

  type <- match.arg(type)
  alternative <- match.arg(alternative)

  if(!is.null(funX))
    stopifnot(is.fv(funX))

  pwrong <- NULL
  use.weights <- !is.null(weights)
  cheat <- !is.null(precomputed)

  if(is.null(argvals) && is.null(observed) && !is.null(funX)) {
    ## assume funX is summary function for observed data
    argvals <- with(funX, .x)
    observed <- with(funX, .y)
    argname.orig <- fvnames(funX, ".x")
    if(is.null(argname)) argname <- argname.orig
    if(is.null(arg.desc))
      arg.desc <- attr(funX, "desc")[match(argname.orig, colnames(funX))]
    theory <- if(use.theory) (theory %orifnull% funX[["theo"]]) else NULL
    if(check) stopifnot(nrow(funX) == nrow(Y)) 
  } else {
    ## construct envelope from raw data
    if(is.null(argname)) argname <- "r"
    if(is.null(arg.desc)) arg.desc <- "distance argument r"
    if(check) {
      ## validate vectors of data
      if(is.null(argvals)) stop("argvals must be supplied")
      if(is.null(observed)) stop("observed must be supplied")
      stopifnot(length(argvals) == nrow(Y))
      stopifnot(length(observed) == length(argvals))
    }
  }

  use.theory <- use.theory && !is.null(theory)
  if(use.theory && check) stopifnot(length(theory) == length(argvals))

  simvals <- Y
  fX <- observed

  if(!is.null(funX)) {
    atr <- attributes(funX)
    fname <- atr$fname
    yexp <- atr$yexp
  } else {
    fname <- "f"
    yexp <- substitute(f(r), list(r=as.name(argname)))
    atr <- list(alim=range(argvals),
                ylab=yexp,
                yexp=yexp,
                fname="f")
  }
  ## for use in making labels
  ## doesn't work at present!
  ## FNAME <- fname.orig %orifnull% fname
  ## TRA   <- transform
  FNAME <- fname
  TRA <- NULL
  
  NAvector <- rep(NA_real_, length(argvals))
  
  if(!cheat) {
    ## ................   standard calculation .....................
    ## validate weights
    if(use.weights && !gaveup) 
      check.nvector(weights, ncol(simvals), 
                    things="simulated functions", naok=TRUE, vname="weights")

    ## determine numbers of columns used
    Ncol <- if(!gaveup) ncol(simvals) else Inf
    if(Ncol < 2)
      stop("Need at least 2 columns of function values")

    ## all columns are used unless 'nsim' or 'jsim' given.
    if(!(is.null(nsim) && is.null(jsim))) {
      if(is.null(jsim)) {
        jsim <- 1:nsim
      } else if(is.null(nsim)) {
        nsim <- length(jsim)
      } else stopifnot(length(jsim) == nsim)
      if(nsim > Ncol)
        stop(paste(nsim, "simulations are not available; only",
                   Ncol, "columns provided"))
    }
    
    ## nsim2 or jsim.mean may be given, and imply dual calculation
    if(!(is.null(nsim2) && is.null(jsim.mean))) {
      if(is.null(jsim.mean)) {
        jsim.mean <- setdiff(seq_len(Ncol), jsim)[1:nsim2]
      } else if(is.null(nsim2)) {
        nsim2 <- length(jsim.mean)
      } else stopifnot(length(jsim.mean) == nsim2)
      if(nsim + nsim2 > Ncol)
        stop(paste(nsim, "+", nsim2, "=", nsim+nsim2, 
                   "simulations are not available; only",
                   Ncol, "columns provided"))
      if(length(intersect(jsim, jsim.mean)))
        warning("Internal warning: Indices in jsim and jsim.mean overlap")
    }
      
    restrict.columns <- !is.null(jsim)
    dual <- !is.null(jsim.mean)

  } else {
    ## ................ precomputed values ..................
    ## validate weights
    if(use.weights) 
      check.nvector(weights, nsim,
                    things="simulations", naok=TRUE, vname="weights")
    restrict.columns <- FALSE
    dual <- FALSE
  }

  shadenames <- NULL
  nsim.mean <- NULL
  
  switch(type,
         pointwise = {
           ## ....... POINTWISE ENVELOPES ...............................
           if(gaveup) {
             lo <- hi <- NAvector
           } else if(cheat) {
             stopifnot(checkfields(precomputed, c("lo", "hi")))
             lo <- precomputed$lo
             hi <- precomputed$hi
           } else {
             simvals[is.infinite(simvals)] <- NA
             if(restrict.columns) {
               simvals <- simvals[,jsim]
               if(use.weights) weights <- weights[jsim]
             }
             nsim <- ncol(simvals)
             if(nrank == 1L) {
               lohi <- apply(simvals, 1L, range)
             } else {
               lohi <- apply(simvals, 1L,
#                             function(x, n) { sort(x)[n] },
                             orderstats,
                             k=c(nrank, nsim-nrank+1L))
             }
             lo <- lohi[1L,]
             hi <- lohi[2L,]
           }
           lo.name <- "lower pointwise envelope of %s from simulations"
           hi.name <- "upper pointwise envelope of %s from simulations"
           ##
           if(!gaveup)
             switch(alternative,
                    two.sided = { },
                    less = {
                      hi <- rep.int(Inf, length(hi))
                      hi.name <- "infinite upper limit"
                    },
                    greater = {
                      lo <- rep.int(-Inf, length(lo))
                      lo.name <- "infinite lower limit"
                    })
           if(use.theory) {
             results <- data.frame(r=argvals,
                                   obs=fX,
                                   theo=theory,
                                   lo=lo,
                                   hi=hi)
           } else {
             m <- if(gaveup) NAvector else
                  if(cheat) precomputed$mmean else 
                  if(!use.weights) apply(simvals, 1L, mean, na.rm=TRUE) else
                  apply(simvals, 1L, weighted.mean, w=weights, na.rm=TRUE)
             results <- data.frame(r=argvals,
                                   obs=fX,
                                   mmean=m,
                                   lo=lo,
                                   hi=hi)
           }
           colnames(results)[1] <- argname
           shadenames <- c("lo", "hi")
           if(do.pwrong) {
             ## estimate the p-value for the 'wrong test'
             if(gaveup) {
               pwrong <- NA_real_
             } else if(cheat) {
               pwrong <- precomputed$pwrong
               do.pwrong <- !is.null(pwrong) && !badprobability(pwrong, FALSE)
             } else {
               dataranks <- t(apply(simvals, 1, rank, ties.method="random"))
               upper.signif <- (dataranks <= nrank)
               lower.signif <- (dataranks >= nsim-nrank+1L)
               is.signif <- switch(alternative,
                                   less = lower.signif,
                                   greater = upper.signif,
                                   two.sided = lower.signif | upper.signif)
               is.signif.somewhere <- matcolany(is.signif)
               pwrong <- sum(is.signif.somewhere)/nsim
             }
           }
         },
         global = {
           ## ..... SIMULTANEOUS ENVELOPES ..........................
           if(gaveup) {
             lo <- hi <- reference <- NAvector
           } else if(cheat) {
             ## ... use precomputed values ..
             stopifnot(checkfields(precomputed, c("lo", "hi")))
             lo <- precomputed$lo
             hi <- precomputed$hi
             if(use.theory) {
               reference <- theory
             } else {
               stopifnot(checkfields(precomputed, "mmean"))
               reference <- precomputed$mmean
             }
             domain <- rep.int(TRUE, length(argvals))
           } else {
             ## ... normal case: compute envelopes from simulations
             if(!is.null(ginterval)) {
               domain <- (argvals >= ginterval[1L]) & (argvals <= ginterval[2L])
               funX <- funX[domain, ]
               simvals <- simvals[domain, ]
             } else domain <- rep.int(TRUE, length(argvals))
             simvals[is.infinite(simvals)] <- NA
             if(use.theory) {
               reference <- theory[domain]
               if(restrict.columns) {
                 simvals <- simvals[, jsim]
                 if(use.weights) weights <- weights[jsim]
               }
             } else if(dual) {
               # Estimate the mean from one set of columns
               # Form envelopes from another set of columns
               simvals.mean <- simvals[, jsim.mean]
               # mmean <-
               reference <- 
                 if(!use.weights) apply(simvals.mean, 1L, mean, na.rm=TRUE) else
                 apply(simvals.mean, 1L, weighted.mean, w=weights[jsim.mean],
                       na.rm=TRUE)
               nsim.mean <- ncol(simvals.mean)
               # retain only columns used for envelope
               simvals <- simvals[, jsim]
             } else {
               # Compute the mean and envelopes using the same data
               if(restrict.columns) {
                 simvals <- simvals[, jsim]
                 if(use.weights) weights <- weights[jsim]
               }
               # mmean <-
               reference <- 
                 if(!use.weights) apply(simvals, 1L, mean, na.rm=TRUE) else
                 apply(simvals, 1L, weighted.mean, w=weights, na.rm=TRUE)
             }
             nsim <- ncol(simvals)
             # compute deviations
             deviations <- sweep(simvals, 1L, reference)
             deviations <-
               switch(alternative,
                      two.sided = abs(deviations),
                      greater = if(clamp) pmax(0, deviations) else deviations,
                      less = if(clamp) pmax(0, -deviations) else (-deviations))
             deviations <- matrix(deviations,
                                  nrow=nrow(simvals), ncol=ncol(simvals))
             ## rescale ?
             sc <- 1
             if(!is.null(scale)) {
               stopifnot(is.function(scale))
               sc <- scale(argvals)
               sname <- paste0("scale", paren(argname))
               tname <- paste("values of", argname)
               ans <- check.nvector(sc, length(argvals), things=tname,
                                    fatal=FALSE, vname=sname)
               if(!ans)
                 stop(attr(ans, "whinge"), call.=FALSE)
               if(any(bad <- (sc <= 0))) {
                 ## issue a warning unless this only happens at r=0
                 if(any(bad[argvals > 0]))
                   warning(paste("Some values of", sname,
                                 "were negative or zero:",
                                 "scale was reset to 1 for these values"),
                           call.=FALSE)
                 sc[bad] <- 1
               }
               deviations <- sweep(deviations, 1L, sc, "/")
             }
             ## compute max (scaled) deviations
             suprema <- apply(deviations, 2L, max, na.rm=TRUE)
             # ranked deviations
             dmax <- sort(suprema)[nsim-nrank+1L]
             # simultaneous bands
             lo <- reference - sc * dmax
             hi <- reference + sc * dmax
           }

           lo.name <- "lower critical boundary for %s"
           hi.name <- "upper critical boundary for %s"

           if(!gaveup)
             switch(alternative,
                    two.sided = { },
                    less = {
                      hi <- rep.int(Inf, length(hi))
                      hi.name <- "infinite upper boundary"
                    },
                    greater = {
                      lo <- rep.int(-Inf, length(lo))
                      lo.name <- "infinite lower boundary"
                    })

           if(use.theory) {
             results <- data.frame(r=argvals[domain],
                                   obs=fX[domain],
                                   theo=reference,
                                   lo=lo,
                                   hi=hi)
           } else {
             results <- data.frame(r=argvals[domain],
                                   obs=fX[domain],
                                   mmean=reference,
                                   lo=lo,
                                   hi=hi)
           }
           colnames(results)[1] <- argname
           shadenames <- c("lo", "hi")
           if(do.pwrong)
             warning(paste("Argument", sQuote("do.pwrong=TRUE"), "ignored;",
                           "it is not relevant to global envelopes"))
         },
         variance={
           ## ....... POINTWISE MEAN, VARIANCE etc ......................
           if(gaveup) {
             Ef <- varf <- NAvector
           } else if(cheat) {
             # .... use precomputed values ....
             stopifnot(checkfields(precomputed, c("Ef", "varf")))
             Ef   <- precomputed$Ef
             varf <- precomputed$varf
           } else {
             ## .... normal case: compute from simulations
             simvals[is.infinite(simvals)] <- NA
             if(restrict.columns) {
               simvals <- simvals[, jsim]
               if(use.weights) weights <- weights[jsim]
             }
             nsim <- ncol(simvals)
             if(!use.weights) {
               Ef   <- apply(simvals, 1L, mean, na.rm=TRUE)
               varf <- apply(simvals, 1L, var,  na.rm=TRUE)
             } else {
               Ef   <- apply(simvals, 1L, weighted.mean, w=weights, na.rm=TRUE)
               varf <- apply(simvals, 1L, weighted.var,  w=weights, na.rm=TRUE)
             }
           }
           if(gaveup) {
             sd <- stdres <- lo <- hi <- loCI <- hiCI <- NAvector
           } else {
             ## derived quantities
             sd <- sqrt(varf)
             stdres <- (fX-Ef)/sd
             stdres[!is.finite(stdres)] <- NA
             ## critical limits
             lo <- Ef - nSD * sd
             hi <- Ef + nSD * sd
             ## confidence interval 
             loCI <- Ef - nSD * sd/sqrt(nsim)
             hiCI <- Ef + nSD * sd/sqrt(nsim)
           }
           lo.name <- paste("lower", nSD, "sigma critical limit for %s")
           hi.name <- paste("upper", nSD, "sigma critical limit for %s")
           loCI.name <- paste("lower", nSD, "sigma confidence bound",
                              "for mean of simulated %s")
           hiCI.name <- paste("upper", nSD, "sigma confidence bound",
                                "for mean of simulated %s")
           ##
           if(!gaveup)
             switch(alternative,
                    two.sided = { },
                    less = {
                      hi <- hiCI <- rep.int(Inf, length(hi))
                      hi.name <- "infinite upper boundary"
                      hiCI.name <- "infinite upper confidence limit"
                    },
                    greater = {
                      lo <- loCI <- rep.int(-Inf, length(lo))
                      lo.name <- "infinite lower boundary"
                      loCI.name <- "infinite lower confidence limit"
                    })
           ## put together
           if(use.theory) {
             results <- data.frame(r=argvals,
                                   obs=fX,
                                   theo=theory,
                                   lo=lo,
                                   hi=hi)
             colnames(results)[1] <- argname
             shadenames <- c("lo", "hi")
             morestuff <- data.frame(mmean=Ef,
                                     var=varf,
                                     res=fX-Ef,
                                     stdres=stdres,
                                     loCI=loCI,
                                     hiCI=hiCI)
             loCIlabel <-
               if(alternative == "greater" && !gaveup) "-infinity" else
               makefvlabel(NULL, NULL, FNAME, "loCI", argname=argname, pre=TRA)
             hiCIlabel <-
               if(alternative == "less" && !gaveup) "infinity" else 
               makefvlabel(NULL, NULL, FNAME, "hiCI", argname=argname, pre=TRA)
             mslabl <- c(
               makefvlabel(NULL, "bar", FNAME, argname=argname, pre=TRA),
               makefvlabel("var", "hat", FNAME, argname=argname, pre=TRA),
               makefvlabel("res", "hat", FNAME, argname=argname, pre=TRA),
               makefvlabel("stdres", "hat", FNAME, argname=argname, pre=TRA),
               loCIlabel,
               hiCIlabel)
             wted <- if(use.weights) "weighted" else NULL
             msdesc <- c(pasteN(wted, "sample mean of %s from simulations"),
                         pasteN(wted, "sample variance of %s from simulations"),
                         "raw residual",
                         "standardised residual",
                         loCI.name, hiCI.name)
           } else {
             results <- data.frame(r=argvals,
                                   obs=fX,
                                   mmean=Ef,
                                   lo=lo,
                                   hi=hi)
             colnames(results)[1] <- argname
             shadenames <- c("lo", "hi")
             morestuff <- data.frame(var=varf,
                                     res=fX-Ef,
                                     stdres=stdres,
                                     loCI=loCI,
                                     hiCI=hiCI)
             loCIlabel <-
               if(alternative == "greater" && !gaveup) "-infinity" else
               makefvlabel(NULL, NULL, FNAME, "loCI", argname=argname, pre=TRA)
             hiCIlabel <-
               if(alternative == "less" && !gaveup) "infinity" else 
               makefvlabel(NULL, NULL, FNAME, "hiCI", argname=argname, pre=TRA)
             mslabl <- c(
               makefvlabel("var", "hat", FNAME, argname=argname, pre=TRA),
               makefvlabel("res", "hat", FNAME, argname=argname, pre=TRA),
               makefvlabel("stdres", "hat", FNAME, argname=argname, pre=TRA),
               loCIlabel,
               hiCIlabel)
             msdesc <- c(pasteN(if(use.weights) "weighted" else NULL,
                                "sample variance of %s from simulations"),
                         "raw residual",
                         "standardised residual",
                         loCI.name, hiCI.name)
           }
           if(do.pwrong) {
             ## estimate the p-value for the 'wrong test'
             if(gaveup) {
               pwrong <- NA_real_
             } else if(cheat) {
               pwrong <- precomputed$pwrong
               do.pwrong <- !is.null(pwrong) && !badprobability(pwrong, FALSE)
             } else {
               upper.signif <- (simvals > hi)
               lower.signif <- (simvals < lo)
               is.signif <- switch(alternative,
                                   less = lower.signif,
                                   greater = upper.signif,
                                   two.sided = lower.signif | upper.signif)
#               is.signif.somewhere <- apply(is.signif, 2, any)
               is.signif.somewhere <- matcolany(is.signif)
               pwrong <- sum(is.signif.somewhere)/nsim
             }
           }
         }
         )

  ############  WRAP UP #########################

  if(use.theory) {
    # reference is computed curve `theo'
    reflabl <- makefvlabel(NULL, NULL, FNAME, "theo", argname=argname, pre=TRA)
    refdesc <- pasteN(theory.adjective,
                      "theoretical value of %s",
                      if(csr) "for CSR" else NULL)
  } else {
    # reference is sample mean of simulations
    reflabl <- makefvlabel(NULL, "bar", FNAME, argname=argname, pre=TRA)
    refdesc <- pasteN(if(use.weights) "weighted" else NULL,
                      "sample mean of %s from simulations")
  }

  lolabl <- if(alternative == "greater" && !gaveup) "-infinity" else
             makefvlabel(NULL, "hat", FNAME, "lo", argname=argname, pre=TRA)
  hilabl <- if(alternative == "less"&& !gaveup) "infinity" else
             makefvlabel(NULL, "hat", FNAME, "hi", argname=argname, pre=TRA)

  result <- fv(results,
               argu=argname,
               ylab=atr$ylab,
               valu="obs",
               fmla= paste(". ~", argname),
               alim=intersect.ranges(atr$alim, range(results[[argname]])),
               labl=c(argname,
                      makefvlabel(NULL, "hat", FNAME,
                                  "obs", argname=argname, pre=TRA),
                      reflabl,
                      lolabl,
                      hilabl),
               desc=c(arg.desc,
                 "observed value of %s for data pattern",
                 refdesc, lo.name, hi.name),
               fname=atr$fname,
               yexp =atr$yexp)

  # columns to be plotted by default
  dotty <- c("obs", if(use.theory) "theo" else "mmean", "hi", "lo")

  if(type == "variance") {
    # add more stuff
    result <- bind.fv(result, morestuff, mslabl, msdesc)
    if(use.theory) dotty <- c(dotty, "mmean")
  }

  fvnames(result, ".") <- dotty
  fvnames(result, ".s") <- shadenames

  unitname(result) <- unitname(funX)
  class(result) <- c("envelope", class(result))

  # tack on envelope information
  attr(result, "einfo") <- list(global = (type =="global"),
                                ginterval = ginterval,
                                alternative=alternative,
                                scale = scale,
                                clamp = clamp,
                                csr = csr,
                                use.theory = use.theory,
                                theory.adjective = theory.adjective,
                                csr.theo = csr && use.theory,
                                simtype = "funs",
                                constraints = "",
                                nrank = nrank,
                                nsim = nsim,
                                VARIANCE = (type == "variance"),
                                nSD = nSD,
                                valname = NULL,
                                dual = dual,
                                nsim = nsim,
                                nsim2 = nsim.mean,
                                Yname = Yname,
                                do.pwrong=do.pwrong,
                                use.weights=use.weights,
                                gaveup = gaveup)

  # tack on saved functions
  if(savefuns && !gaveup) {
    nSim <- ncol(Y)
    alldata <- cbind(argvals, Y)
    simnames <- paste("sim", 1:nSim, sep="")
    colnames(alldata) <- c(argname, simnames)
    alldata <- as.data.frame(alldata)
    SimFuns <- fv(alldata,
                  argu=argname,
                  ylab=atr$ylab,
                  valu="sim1",
                  fmla= paste(". ~", argname),
                  alim=atr$alim,
                  labl=names(alldata),
                  desc=c(arg.desc,
                          paste("Simulation ", 1:nSim, sep="")),
                  unitname=unitname(funX))
    fvnames(SimFuns, ".") <- simnames
    attr(result, "simfuns") <- SimFuns
  }
  if(do.pwrong)
    attr(result, "pwrong") <- pwrong
  if(use.weights)
    attr(result, "weights") <- weights
  return(result)
}
