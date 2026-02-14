#'
#' bw.pcf.R
#'
#' $Revision: 1.11 $  $Date: 2026/02/14 10:06:42 $
#'
#' bandwidth selection for pcf
#' with least-squares cross-validation method
#' 
#' Original code by: Rasmus Waagepetersen and Abdollah Jalilian
#'
#' Extended version by: Adrian Baddeley, Tilman Davies and Martin Hazelton
#' 
#' Guan, Y. (2007). A composite likelihood cross-validation approach in 
#'   selecting bandwidth for the estimation of the pair correlation function. 
#'   Scandinavian Journal of Statistics, 34(2), 336--346. 
#'   DOI: http://doi.org/10.1111/j.1467-9469.2006.00533.x
#' Guan, Y. (2007). A least-squares cross-validation bandwidth 
#'   selection approach in pair correlation function estimations. 
#'   Statistics & Probability Letters, 77(18), 1722--1729. 
#'   DOI: http://doi.org/10.1016/j.spl.2007.04.016
#' Jalilian, A. and Waagepetersen, R. (2018) 
#'   Fast bandwidth selection for estimation of the pair correlation function.
#'   Journal of Statistical Computation and Simulation 88(10) 2001--2011.
#'   DOI: 10.1080/00949655.2018.1428606
#' Baddeley, A, Davies, T.M. and Hazelton, M.L. (2025)
#'   An improved estimator of the pair correlation function
#'   of a spatial point process. Biometrika, to appear.
#'
#' Copyright (c) 2017-2025 Rasmus Waagepetersen, Abdollah Jalilian,
#'               Adrian Baddeley, Tilman Davies and Martin Hazelton

bw.pcf <- function(X, ..., rmax=NULL, nr=10000, 
                      cv.method=c("compLik", "leastSQ", "oracle"),
                      leaveoneout=TRUE, simple=TRUE,
                      fast=TRUE, srange=NULL, ns=32,
                      use.count=TRUE, gtrue=NULL, 
                      verbose=FALSE, warn=TRUE)
{
  cv.method <- match.arg(cv.method)
  BWPCFengine(X, ...,
                 rmax=rmax, nr=nr, cv.method=cv.method,
                 leaveoneout=leaveoneout, simple=simple,
                 fast=fast, srange=srange, ns=ns,
                 use.count=use.count,
                 gtrue=gtrue, verbose=verbose, warn=warn)
}

bw.pcfinhom <-
  function(X, lambda=NULL, ..., rmax=NULL, nr=10000, 
           cv.method=c("compLik", "leastSQ", "oracle"),
           leaveoneout=TRUE, simple=TRUE,
           fast=TRUE, srange=NULL, ns=32,
           use.count=TRUE, gtrue=NULL, 
           verbose=FALSE, warn=TRUE)
{
  cv.method <- match.arg(cv.method)
  BWPCFengine(X, ..., lambda=lambda,
                 rmax=rmax, nr=nr, cv.method=cv.method,
                 leaveoneout=leaveoneout, simple=simple,
                 fast=fast, srange=srange, ns=ns,
                 use.count=use.count,
                 gtrue=gtrue, verbose=verbose, warn=warn)
}

BWPCFengine <- function(X, ..., lambda=NULL,
                           rmax=NULL, nr=10000,
                           correction=c("translate", "isotropic"),
                           divisor=c("r", "d", "a", "t"),
                           kernel = "epanechnikov",
                           zerocor=c("none", "weighted", "convolution",
                                     "reflection", "bdrykern", "JonesFoster"),
                           gref=NULL,
                           cv.method=c("compLik", "leastSQ", "martin", "oracle"),
                           gtrue=NULL, 
                           adaptive=FALSE,
                           leaveoneout=TRUE, simple=TRUE,
                           fast=TRUE, srange=NULL, ns=32,
                           use.count=TRUE,
                           renormalise=TRUE, normpower=1,
                           verbose=FALSE, debug=FALSE,
                           warn=TRUE)
{
  stopifnot(is.ppp(X))
  X <- unmark(X)
  win <- Window(X)
  areaW <- area(win)
  nX <- npoints(X)

  if(cv.method == "oracle") stopifnot(is.function(gtrue))

  homogeneous <- is.null(lambda)
  
  if(adaptive) warning("adaptive=TRUE is currently ignored!")

  if(missing(correction)) {
    correction <- match.arg(correction)
  } else {
    correction <- pickoption("correction", correction,
                             c(isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               good="translate",
                               best="isotropic"),
                             multi=FALSE)
  }

  divisor <- match.arg(divisor)
  kernel <- match.kernel(kernel)
  zerocor <- match.arg(zerocor)
  cv.method <- match.arg(cv.method)

  #' maximum distance lag: rmax
  rmax.max     <- diameter(Frame(X))
  rmax.default <- rmax.rule("K", win,  nX/areaW)
  if(is.null(rmax)) {
    rmax <- if(cv.method == "martin") rmax.max else rmax.default
  } else {
    check.1.real(rmax)
    rmax <- min(rmax, rmax.max)
  }
  #' number of subintervals for discretization of [0, rmax]: nr
  #' length of subintervals
  discr <- rmax / nr
  #' breaks of subintervals
  rs <- seq(0, rmax, length.out= nr + 1)

  #' range of bandwidths to try
  if(is.null(srange))
    srange <- c(0, rmax.default/4)
  
  #' closepairs distances: \\ u - v \\
  #' Pre-compute close pair distances for use in 'pcf'
  #'   we need close pairs up to a distance rmax + smax
  #'   where 'smax' is the maximum halfwidth of the support of the kernel
  smax <- srange[2] * (if(kernel == "gaussian") 2 else kernel.factor(kernel))
  cpfull <- closepairs(X, rmax + smax, what="all", twice=TRUE)
  
  #' For cross-validation, restrict close pairs to distance rmax 
  ok <- (cpfull$d <= rmax)
  cp <- lapply(cpfull, "[", i=ok)

  ds <- cp$d
  npairs <- length(ds)
  #' determining closepairs distances are in which subinterval
  idx <- round(ds / discr) + 1L
  idx <- pmin.int(idx, nr+1L)
  
  #' edge correction factor associated with each (i,j)
  edgewt <-
    if(leaveoneout) {
      switch(correction,
             translate = edge.Trans(dx=cp$dx, dy=cp$dy, W=win, paired=TRUE),
             isotropic = edge.Ripley(X[cp$i], ds, W=win))
    } else 1

  #' normalisation factors
  if(homogeneous) {
    #' homogeneous case
    lambdaBar <- nX/areaW
    lambda2area <- lambdaBar^2 * areaW
    renorm.factor <- 1
    #' expected value of \sum 1/g(d_{ij})
    expectedsum <- if(cv.method != "martin") NA else
                   if(use.count) npairs else
                   (lambdaBar^2 * integral(setcov(win), domain=disc(rmax)))
    #'
    inhomargs <- list()
  } else {
    #' inhomogeneous case
    a <- resolve.reciplambda(X, lambda=lambda, ..., leaveoneout=leaveoneout)
    lambdaX      <- a$lambda
    reciplambdaX <- a$reciplambda
    edgewt <- edgewt * reciplambdaX[cp$i] * reciplambdaX[cp$j]
    lambda2area <- areaW
    if(renormalise && nX > 0) {
      check.1.real(normpower)
      stopifnot(normpower %in% 1:2)
      renorm.factor <- (areaW/sum(reciplambdaX))^normpower
    } else {
      renorm.factor <- 1
    }
    expectedsum <- if(cv.method != "martin") NA else
                   stop("cv='martin' is not supported for bw.pcfinhom")
    inhomargs <- list(lambda=lambdaX,
                      renormalise=renormalise,
                      normpower=normpower)
  } 

  #' arguments to pcf or pcfinhom
  pcfargs <- resolve.defaults(
    list(X=X,
         r=rs,
         kernel=kernel,
         correction=correction,
         divisor=divisor,
         zerocor=zerocor,
         gref = gref,
         close=cpfull),
    list(...),
    inhomargs)

  ## evaluate true pcf at 'rs' if known
  gtruers <- if(cv.method == "oracle") gtrue(rs) else NULL

  stuff <- list(cv.method=cv.method,
                homogeneous=homogeneous,
                nX=nX,
                rs=rs,
                cp=cp,
                ds=ds,
                idx=idx,
                edgewt=edgewt,
                simple=simple,
                expectedsum=expectedsum,
                pcfargs=pcfargs,
                lambda=lambda,
                lambda2area=lambda2area,
                renorm.factor=renorm.factor,
                discr=discr,
                gtrue=gtrue,
                gtruers=gtruers,
		show=verbose || debug,
                debug=debug)
  stuff <- list2env(stuff)

  #' find optimum bandwidth
  optimum <- switch(cv.method,
                    compLik = "max",
                    leastSQ  = "min",
                    martin  = "min",
                    oracle = "min")
  if(fast) {
    #' optimize
    z <- optimizeWithTrace(CVforPCF, srange, stuff=stuff,
                           maximum=(optimum == "max"))
    ox <- order(z$x)
    sigma  <- z$x[ox]
    cv     <- z$y[ox]
  } else {
    #' evaluate at an evenly-spaced grid of values
    check.1.integer(ns)
    stopifnot(ns > 2)
    sigma <- seq(srange[1], srange[2], length.out=ns)
    cv <- numeric(ns)
    valid <- (sigma > 0)
    cv[!valid] <- if(optimum == "max") -Inf else Inf
    for(i in which(valid)) {
      cv[i] <- CVforPCF(sigma[i], stuff)
    }
  }
  #' pack up
  criterion <- switch(cv.method,
                      compLik = "composite likelihood cross-validation",
                      leastSQ = "least squares cross-validation",
                      martin = "Hazelton")
  result <- bw.optim(cv, sigma, optimum=optimum,
                     criterion = criterion,
                     warnextreme=warn, hargnames=c("rmax", "srange"),
                     unitname=unitname(X))
  return(result)
}

.Integrand <- function(x, g) { 2 * pi * x * g(x) }

CVforPCF <- function(BW, stuff) {
  #' Compute cross-validation objective function
  force(BW)
  #' (1) Form estimates of pcf using bandwidth 'BW'
  if(stuff$pcfargs$zerocor != "JonesFoster") {
    #' usual case
    A <- PCFatRandD(stuff, bw=BW)
    stuff$grs <- A$grs  # estimates g(r) at grid of 'r' values
    stuff$gds <- A$gds  # estimates g(d) at pairwise distances 'd'
  } else {
    #' Jones-Foster estimator
    Aconv <- PCFatRandD(stuff, bw=BW, zerocor="convolution")
    Abdry <- PCFatRandD(stuff, bw=BW, zerocor="bdrykern")
    stuff$grs <- Aconv$grs * exp(Abdry$grs/Aconv$grs - 1)
    stuff$gds <- Aconv$gds * exp(Abdry$gds/Aconv$gds - 1)
  }
  with(stuff, {
    #' (2) remove negative and zero values and NaN etc
    if(debug) gds.save <- gds
    grs[!is.finite(grs)] <- 0
    gds[!is.finite(gds)] <- .Machine$double.eps
    grs <- pmax.int(grs, 0)
    gds <- pmax.int(gds, .Machine$double.eps)
    if(debug) {
      plot(gds.save, gds, xlab="g(d_[ij])", ylab="g(d_[ij])^{-ij}")
      abline(0,1)
      splat("range g_after/g_before = ", prange(range(gds/gds.save)))
    }
    #' (3) compute value of objective function
    switch(cv.method,
           compLik={
             #' composite likelihood cross-validation
             #' the integral term: 2 \pi \int_{0}^{rmax} \hat g(r) r dr
             normconst <- 2 * pi * sum(grs * rs) * discr
             value <- mean(log(gds)) - log(normconst)
             if(debug)
               splat("normconst = ", normconst)
           },
           leastSQ={
             #' least squares cross-validation
             #' the integral term: 2 \pi \int_{0}^{rmax} \hat g^2(r) r dr
             normconst <- 2 * pi * sum(grs^2 * rs) * discr
             value <- normconst - 2 * sum(gds * edgewt / (lambda2area))
             if(debug)
               splat("normconst = ", normconst)
           },
           martin = {
             #' Martin Hazelton's criterion
             value <- (sum(1/gds) - expectedsum)^2
           },
           oracle = {
             value <- discr * sum((grs - gtruers)^2)
           },
           stop("Unrecognised cross-validation method"))
    #' debug
    if(debug && !is.finite(value)) {
      cat("grs:")
      print(summary(gds))
      cat("gds:")
      print(summary(gds))
      splat("discr:", discr,
            "sum(grs * rs):", sum(grs*rs),
            "normconst:", normconst)
    }
    if(show || !is.finite(value))
      splat("Returning cv", paren(BW, "["), "=", value)
    return(value)
  })
}

PCFatRandD <- function(stuff, bw=NULL, zerocor=NULL) {
  #' compute estimates of pcf
  #'     grs = g(rs) at a grid of distances 'rs'
  #' and
  #      gds = g(ds) at the pairwise distances 'ds'
  #' optionally applying leave-one-out rule.
  #'
  #' Optionally override settings in 'stuff' 
  stuff$BW <- bw %orifnull% stuff$bw
  stuff$ZEROCOR <- zerocor %orifnull% stuff$pcfargs$zerocor
  stuff$KERNEL <- stuff$pcfargs$kernel
  #'
  #' Compute ...
  with(stuff, {
    if(show) splat("BW=", BW)
    #' values of pair correlation at breaks of subintervals
    a <- append(pcfargs, list(bw=BW))
    gfun <- if(homogeneous) do.call(pcf, a) else do.call(pcfinhom, a)
    if(show) plot(gfun)
    grs <- gfun[[fvnames(gfun, ".y")]]
    #' make sure that the estimated pair correlation at origin is finite
    if (!is.finite(grs[1]))
      grs[1] <- grs[2]
    #' approximate the pair correlation values at closepairs distances
    gds <- grs[idx]
    if(debug) gds.save <- gds
    if(leaveoneout) {
      #' Remove pairs to approximate the cross-validation term: g^{-(u, v)}
      switch(divisor,
             r = {
               #' Subtract contribution 
               if (simple) {
                 wt <- (edgewt / (2 * pi * ds * lambda2area)) * renorm.factor
                 gds <- gds - 2 * wt * dkernelBC(ds, ds,
                                                 sd=BW, kernel=KERNEL,
                                                 zerocor=ZEROCOR)
               } else {
                 cpi <- cp$i
                 cpj <- cp$j
                 ww <- (edgewt / (2 * pi * lambda2area)) * renorm.factor
                 for (k in 1:length(ds)) {
                   #' adjust estimate at r = ds[k]
                   exclude <- (cpi == cpi[k]) | (cpj == cpj[k])
                   kernex <- dkernelBC(ds[k], ds[exclude], 
                                       sd=BW, kernel=KERNEL, zerocor=ZEROCOR)
                   gds[k] <- gds[k] - 2 * sum(ww[exclude] * kernex)/ds[k] 
                 }
               }
             },
             d = {
               #' Compute weight associated with each (i,j) pair
               #' Note divisor='r' and divisor='d' have the same effect at r=d
               wt <- (edgewt / (2 * pi * ds * lambda2area)) * renorm.factor
               #' Subtract contribution 
               if (simple) {
                 gds <- gds - 2 * wt * dkernelBC(ds, ds,
                                                 sd=BW, kernel=KERNEL,
                                                 zerocor=ZEROCOR)
               } else {
                 cpi <- cp$i
                 cpj <- cp$j
                 for (k in 1:length(ds)) {
                   exclude <- (cpi == cpi[k]) | (cpj == cpj[k])
                   kernex <- dkernelBC(ds[k], ds[exclude], 
                                       sd=BW, kernel=KERNEL, zerocor=ZEROCOR)
                   gds[k] <- gds[k] - 2 * sum(wt[exclude] * kernex)
                 }
               }
             },
             a = ,
             t = {
               #' Compute weight associated with each (i,j) pair
               wt <- (edgewt / lambda2area) * renorm.factor
               #' convert to transformed scale
               if(divisor == "a") {
                 #' convert to areas
                 As <- pi * ds^2
               } else {
                 #' transformation determined by integral 2 pi x g(x) dx
                 As <- indefinteg(.Integrand, ds, g=gref, lower=0)
               }
               #' extract bandwidth used on area scale
               Abw <- attr(gfun, "bw")
               #' Subtract contribution 
               if (simple) {
                 gds <- gds - 2 * wt * dkernelBC(As, As,
                                                 sd=Abw, kernel=KERNEL,
                                                 zerocor=ZEROCOR)
               } else {
                 cpi <- cp$i
                 cpj <- cp$j
                 for (k in 1:length(ds)) {
                   exclude <- (cpi == cpi[k]) | (cpj == cpj[k])
                   kernex <- dkernelBC(As[k], As[exclude],
                                       sd=Abw, kernel=KERNEL,
                                       zerocor=ZEROCOR)
                   gds[k] <- gds[k] - 2 * sum(wt[exclude] * kernex)
                 }
               }
             })
    }
    #' remove negative and zero values
    gds <- pmax.int(.Machine$double.eps, gds)
    if(debug) {
      plot(gds.save, gds, xlab="g(d_[ij])", ylab="g(d_[ij])^{-ij}")
      abline(0,1)
      splat("range g_after/g_before = ", prange(range(gds/gds.save)))
      cat("edgewt:\n")
      print(summary(edgewt))
    }
    #' return
    return(list(grs=grs, gds=gds))
  })
}

