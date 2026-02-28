#
#   envelope.R
#
#   computes simulation envelopes 
#
#   $Revision: 2.133 $  $Date: 2026/02/28 05:12:17 $
#


envelope <- function(Y, fun, ...) {
  UseMethod("envelope")
}

  # .................................................................
  #     A 'simulation recipe' contains the following variables
  #
  #  type = Type of simulation
  #           "csr": uniform Poisson process
  #           "rmh": simulated realisation of fitted Gibbs or Poisson model 
  #          "kppm": simulated realisation of fitted cluster model 
  #          "dppm": simulated realisation of fitted determinantal model 
  #          "slrm": simulated realisation of fitted spatial logistic regression
  #       "process": simulated realisation of given point process 
  #          "expr": result of evaluating a user-supplied expression
  #          "list": user-supplied list of point patterns
  #
  #  expr = expression that is repeatedly evaluated to generate simulations
  #
  #    envir = environment in which to evaluate the expression `expr'
  #
  #    'csr' = TRUE iff the model is (known to be) uniform Poisson
  #
  #    pois  = TRUE if model is known to be Poisson
  #  
  #  constraints = additional information about simulation (e.g. 'with fixed n')
  #
  # ...................................................................

simulrecipe <- function(type, expr, envir, csr, pois=csr, constraints="") {
  if(csr && !pois) warning("Internal error: csr=TRUE but pois=FALSE")
  out <- list(type=type,
              expr=expr,
              envir=envir,
              csr=csr,
              pois=pois,
              constraints=constraints)
  class(out) <- "simulrecipe"
  out
}

make.simulrecipe <- function(object, envir, ...) {
  UseMethod("make.simulrecipe")
}

## //////////////// METHODS FOR "ppp" //////////////////////////
                                 
envelope.ppp <-
  function(Y, fun=Kest, nsim=99, nrank=1, ...,
           funargs=list(), funYargs=funargs,
           simulate=NULL, fix.n=FALSE, fix.marks=FALSE,
           verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE, 
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL,
           maxnerr=nsim, rejectNA=FALSE, silent=FALSE,
           do.pwrong=FALSE,
           envir.simul=NULL) {
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())

  fix.marks <- fix.marks && is.marked(Y)

  #' Data pattern is argument Y
  X <- Y

  if(!is.null(simulate)) {
    #' Simulations are determined by 'simulate' argument
    if(fix.n || fix.marks) 
      warning("fix.n and fix.marks were ignored, because 'simulate' was given")
    #' Processing is deferred to envelopeEngine
    simrecipe <- simulate
  } else {
    ## Default: CSR or binomial process
    simrecipe <- make.simulrecipe(X,
                                  envir.here,
                                  fix.n=fix.n, fix.marks=fix.marks)
  }
  
  envelopeEngine(X=X, fun=fun, simul=simrecipe,
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=clipdata,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname,
                 maxnerr=maxnerr, rejectNA=rejectNA, silent=silent,
                 cl=cl, envir.user=envir.user, do.pwrong=do.pwrong)
}

make.simulrecipe.ppp <- function(object, envir, ...,
                                 fix.n=FALSE, fix.marks=FALSE) {
  X <- object

  nX <- npoints(X)
  Xwin <- Window(X)
  Xmarx <- marks(X)
  Xintens <- nX/area(Xwin)
  assign("nX", nX, envir=envir)
  assign("Xwin",    Xwin,    envir=envir)
  assign("Xmarx", Xmarx, envir=envir)
  assign("Xintens", Xintens, envir=envir)

  if(!fix.n && !fix.marks) {
    #' Realisations of complete spatial randomness with lambda = intensity(X)
    simexpr <- if(is.null(Xmarx)) {
      #' unmarked point pattern
      expression(rpoispp(Xintens, win=Xwin))
    } else if(is.null(dim(Xmarx))) {
      #' single column of marks
      expression({
        A <- rpoispp(Xintens, win=Xwin);
        j <- sample(nX, npoints(A), replace=TRUE);
        A %mark% Xmarx[j]
      })
    } else {
      #' multiple columns of marks
      expression({
        A <- rpoispp(Xintens, win=Xwin);
        j <- sample(nX, npoints(A), replace=TRUE);
        A %mark% Xmarx[j, , drop=FALSE]
      })
    }
    # evaluate in 'envir'
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir,
                             csr   = TRUE,
                             pois  = TRUE)
  } else if(fix.marks) {
    # ...................................................
      # Realisations of binomial process
    # with fixed number of points and fixed marks
    # will be generated by runifpoint
    simexpr <- expression(runifpoint(nX, Xwin) %mark% Xmarx)
    # simulation constraints (explanatory string)
    constraints <-
      if(is.multitype(X)) "with fixed number of points of each type" else
                          "with fixed number of points and fixed marks"
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir,
                             csr   = TRUE,
                             pois  = TRUE,
                             constraints = constraints)
  } else {
    # ...................................................
    # Realisations of binomial process
    # will be generated by runifpoint
    simexpr <- if(is.null(Xmarx)) {
      ## unmarked
      expression(runifpoint(nX, Xwin))
    } else if(is.null(dim(Xmarx))) {
      ## single column of marks
      expression({
        A <- runifpoint(nX, Xwin);
        j <- sample(nX, npoints(A), replace=TRUE);
        A %mark% Xmarx[j]
      })
    } else {
      ## multiple columns of marks
      expression({
        A <- runifpoint(nX, Xwin);
        j <- sample(nX, npoints(A), replace=TRUE);
        A %mark% Xmarx[j, ,drop=FALSE]
      })
    }
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir,
                             csr   = TRUE,
                             pois  = TRUE,
                             constraints = "with fixed number of points")
  }
  return(simrecipe)
}

make.simulrecipe.clusterprocess <- function(object, envir, ..., W) {
  #' Generate realisations of cluster process with given parameters
  assign("simprocess", object,     envir=envir)
  assign("simwin",     as.owin(W), envir=envir)
  simexpr <- expression(simulate(simprocess, win=simwin, drop=TRUE))
  simrecipe <- simulrecipe(type = "process",
                           expr = simexpr,
                           envir = envir,
                           csr   = FALSE,
                           pois  = FALSE)
  return(simrecipe)
}



## Code for envelope.ppm, envelope.kppm, envelope.slrm
## is moved to spatstat.model

