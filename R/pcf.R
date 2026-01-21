#'
#'   pcf.R 
#'
#' Calculate pair correlation function from point pattern (pcf.ppp)
#' 
#' $Revision: 1.79 $ $Date: 2026/01/21 06:26:39 $
#'
#' Copyright (c) 2008-2025 Adrian Baddeley, Tilman Davies and Martin Hazelton


pcf <- function(X, ...) {
  UseMethod("pcf")
}

pcf.NAobject <- function(X, ...) { return(NAobject("fv")) }

pcf.ppp <- function(X, ..., r=NULL, rmax=NULL,
                   adaptive=FALSE,
                   kernel="epanechnikov", bw=NULL, h=NULL,
                   bw.args=list(),
                   stoyan=0.15,
                   adjust = 1,
                   correction=c("translate", "Ripley"),
                   divisor=c("r", "d", "a", "t"),
                   zerocor=c("weighted", "reflection", "convolution",
                             "bdrykern", "JonesFoster", "none"),
                   gref=NULL,
                   tau = 0,
                   fast=TRUE,
                   var.approx=FALSE,
                   domain=NULL, ratio=FALSE,
                   close=NULL)
{
  kernel <- match.kernel(kernel)
  if(is.function(divisor)) divisor <- divisor(X)
  divisor <- match.arg(divisor)
  zerocor <- match.arg(zerocor)
  check.1.real(adjust)

  DEBUG <- isTRUE(getOption("debug.smoothpcf"))
  if(DEBUG) {
    splat("pcf")
    started <- proc.time()
  }
  
  ## ...... get point pattern information .......
  verifyclass(X, "ppp")
  win <- Window(X)
  areaW <- area(win)
  npts <- npoints(X)
  samplesize <- npairs <- npts * (npts - 1)
  lambda <- npts/areaW
  lambda2area <- npairs/areaW
  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambda)        

  ## ....... handle argument 'domain' .......................
  if(!is.null(domain)) {
    if(!(divisor %in% c("r", "d")))
      stop("Sorry, option divisor =", sQuote(divisor),
           "is not yet available when 'domain' is given",
           call.=FALSE)
    if(zerocor != 'none')
      stop(paste0("Sorry, option zerocor=", sQuote(zerocor),
                  "is not yet available when 'domain' is given"),
           call.=FALSE)
    ## estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, win))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    # trick pcfdot() into doing it
    indom <- inside.owin(X$x, X$y, domain)
    marx <- factor(indom, levels=c(FALSE,TRUE))
    g <- pcfdot(X %mark% marx,
                i="TRUE",
                r=r,
                correction=correction, kernel=kernel, bw=bw, stoyan=stoyan,
                divisor=divisor,
                ...)
    if(!ratio) {
      ## relabel
      g <- rebadge.fv(g, quote(g(r)), "g")
    } else {
      ## construct ratfv object
      ninside <- sum(indom)
      samplesize <- ninside * (npts-1)
      g <- ratfv(as.data.frame(g), NULL, samplesize,
                 "r", quote(g(r)),
                 "theo", NULL, c(0, rmaxdefault), 
                 attr(g, "labl"), attr(g, "desc"), fname="g",
                 ratio=TRUE)
    }
    unitname(g) <- unitname(X)
    if(var.approx)
      warning("var.approx is not implemented when 'domain' is given")
    return(g)
  }

  ## ......... edge correction .........................
  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             good="translate",
                             best="best",
                             none="none"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, win$type, correction.given)

  ## .... determine bandwidth ......................

  info <- list(kernel=kernel, divisor=divisor, zerocor=zerocor,
               h.given=h, bw.given=bw, adjust=adjust)
  how <- rule <- NULL
  
  if(!is.null(bw) && !is.null(h))
    stop("Arguments 'h' and 'bw' are incompatible", call.=FALSE)

  ## how the bandwidth will be determined
  if(is.null(bw) && is.null(h)) {
    ## default rule
    how <- "rule"
    rule <- bw <- if(!adaptive) "bw.stoyan" else "bw.abram"
  } else if(!is.null(h)) {
    ## h is given
    if(is.character(h)) {
      how <- "rule"
      rule <- tolower(h)
    } else if(is.numeric(h)) {
      stopifnot(length(h) == 1 && h > 0)
      how <- "h"
    } else stop("h should be a numeric value, or a character string")
  } else {
    ## bw is given
    if(is.numeric(bw) && length(bw) == 1) {
      if(bw <= 0) stop("The bandwidth bw should be positive", call.=FALSE)
      how <- "bw"
    } else if(is.character(bw)) {
      how <- "rule"
      rule <- tolower(bw)
    } else if(is.function(bw)) {
      how <- "fun"
      rule <- bw
    } else stop("bw should be a numeric value, a string, or a function",
                call.=FALSE)
  }
  info <- append(info, list(how=how, rule=rule))

  ## bandwidth arguments bw.args may depend on X
  if(is.function(bw.args)) {
    bw.args <- bw.args(X)
    if(!is.list(bw.args))
      stop("When bw.args is a function, it should return a named list",
           call.=FALSE)
  }

  ## now actually determine the bandwidth from X (unless adaptive = TRUE)
  cker <- kernel.factor(kernel)
  switch(how,
         bw = {
           ## bandwidth is determined by numeric value 'bw'
           h <- bw * cker
         },
         h = {
           ## bandwidth is determined by numeric value 'h'
           bw <- h / cker
         },
         fun = {
           if(!adaptive) {
             ## bandwidth selection *function* applied to X
             bwformals <- names(formals(bw))
             xtra <- list(kernel=kernel,
                          correction=correction[1L],
                          divisor=divisor,
                          zerocor=zerocor,
                          adaptive=adaptive,
                          close=close)
             if(!("..." %in% bwformals)) 
               xtra <- xtra[intersect(names(xtra), bwformals)]
             bw.args <- resolve.defaults(bw.args, xtra)
             bw <- do.call(bw, append(list(quote(X)), bw.args))
             bw.args <- list()
             h <- cker * bw
           }
         },
         rule = {
           ## Character argument 'rule' specifies a bandwidth selection rule
           ## handle the spatial statistics rules now
           switch(rule,
                  bw.stoyan = ,
                  stoyan = {
                    ## Stoyan & Stoyan 1995, eq (15.16), page 285
                    ## for Epanechnikov kernel
                    bw <- stoyan/sqrt(5 * lambda)
                    h <- bw * cker
                  },
                  bw.fiksel = ,
                  fiksel = {
                    ## Fiksel (1988)
                    bw <- 0.1/sqrt(lambda)
                    h <- bw * cker
                  })
           ## (bandwidth may still be 'character')
         })

  #' bandwidth may still be 'character' or 'function'
  
  if(is.numeric(bw)) {
    ## absorb the 'adjust' factor now
    bw <- adjust * bw
    h <- adjust * h
    adjust <- 1
    info <- append(info, list(bw.calc=bw, h.calc=h))
  }

  ########## r values ############################
  # handle arguments r and breaks 

  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # arguments for 'density.default' or 'densityAdaptiveKernel.default'
  denargs <- resolve.defaults(list(kernel=kernel, bw=bw, adjust=adjust),
                              bw.args,
                              list(...),
                              list(n=length(r), from=0, to=rmax),
                              .StripNull = TRUE)

  ############### transformation of distances ################

  switch(divisor,
         r = ,
         d = ,
         a = {
           if(!is.null(gref)) {
             warning(paste("Argument gref is ignored when divisor =",
                           dQuote(divisor)), call.=FALSE)
             gref <- NULL
           }
         },
         t = {
           if(is.null(gref)) {
             ## default: the pcf of the Poisson process
             gref <- function(x) { rep.int(1, length(x)) }
           } else {
             ## normal case: user specified reference function or model
             if(inherits(gref, c("kppm", "dppm", "ppm", "slrm",
                                 "detpointprocfamily", "zclustermodel"))) {
               model <- gref
               if(!requireNamespace("spatstat.model")) 
                 stop("The package spatstat.model is required when",
                      "'gref' is a fitted model",
                      call.=FALSE)
               gref <- spatstat.model::pcfmodel(model)
               if(!is.function(gref))
                 stop("Internal error: pcfodel() did not yield a function",
                      call.=FALSE)
             } else if(!is.function(gref)) {
               stop(paste("Argument", sQuote("gref"),
                          "should be a function or a point process model"),
                    call.=FALSE)
             }
           }
           integrand <- function(x, g) { 2 * pi * x * g(x) }
         })

  #################################################
  ## determine an upper bound on pairwise distances that need to be collected
  hmax <- if(is.numeric(h)) h else (2*cker*stoyan/sqrt(lambda))
  sker <- if(kernel == "gaussian") 4 else 1
  dmax <- rmax + sker * hmax
  if(is.numeric(denargs$bw)) {
    ## compute the bandwidth on the transformed scale now
    switch(divisor,
           r = {},
           d = {},
           a = {
             ## convert to bandwidth(s) for areas
             ## (using same rule as in 'sewpcf')
             bw.area <- with(denargs, pi * (from + to) * bw)
             hmax.area <- cker * max(bw.area)
             ## determine the maximum value of area that needs to be observed
             area.max <- pi * rmax^2
             area.max <- area.max + sker * hmax.area
             ## convert back to a bound on distance
             dmax <- max(dmax, sqrt(area.max/pi))
           },
           t = {
             ## use transformation
             midpoint <- with(denargs, (from + to)/2)
             ## compute derivative of transformation at midpoint of interval
             midslope <- 2 * pi * midpoint * gref(midpoint)
             ## convert bandwidth to transformed scale
             bw.trans <- midslope * denargs$bw
             hmax.trans <- cker * max(bw.trans)
             ## Taylor approx to T^{-1}(T(dmax) + hmax)
             dmax.trans <- dmax + sker * hmax.trans/(2 * pi * dmax * gref(dmax))
             dmax <- max(dmax, dmax.trans)
           })
  }
  info <- append(info, list(rmax=rmax, hmax=hmax, dmax=dmax))

  ## Precompute transform ## and inverse transform
  if(!is.null(gref)) {
    rr <- seq(0, dmax, length.out=16384)
    tt <- indefinteg(integrand, rr, g=gref, lower=0)
    Transform <- approxfun(rr, tt, rule=2)
    ## InvTransform <- approxfun(tt, rr, rule=2)
  }
  
  #######################################################
  ## compute pairwise distances up to 'dmax'
  if(npts > 1) {
    needall <- any(correction %in% c("translate", "isotropic"))
    if(is.null(close)) {
      what <- if(needall) "all" else "ijd"
      close <- closepairs(X, dmax, what=what)
    } else {
      #' check 'close' has correct format
      needed <- if(!needall) c("i", "j", "d") else
                 c("i", "j", "xi", "yi", "xj", "yj", "dx", "dy", "d")
      if(any(is.na(match(needed, names(close)))))
        stop(paste("Argument", sQuote("close"),
                   "should have components named",
                   commasep(sQuote(needed))),
             call.=FALSE)
    }
    dIJ <- close$d
  } else {
    undefined <- rep(NaN, length(r))
  }

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  out <- ratfv(df,
               NULL, samplesize,
               "r", quote(g(r)),
               "theo", NULL,
               alim,
               c("r","%s[Pois](r)"),
               c("distance argument r", "theoretical Poisson %s"),
               fname="g",
               ratio=ratio)

  ###### compute #######

  if(DEBUG) {
    elapsed <- proc.time() - started
    splat("pcf: ready to compute estimates after", codetime(elapsed))
  }
  
  bw.used <- bwvalues.used <- NULL
  
  if(any(correction=="none")) {
    #' uncorrected
    if(npts > 1) {
      kdenN <- sewpcf(d=dIJ, w=1,
                         denargs=denargs,
                         lambda2area=lambda2area,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      if(DEBUG) {
        elapsed <- proc.time() - started
        splat("pcf: returned from sewpcf after", codetime(elapsed))
      }
      gN <- kdenN$g
      bw.used <- attr(kdenN, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenN, "bwvalues")
    } else gN <- undefined
    if(DEBUG) {
      elapsed <- proc.time() - started
      splat("Ready to bind.fv after", codetime(elapsed))
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(un=gN),
                      denominator=samplesize,
                      labl="hat(%s)[un](r)",
                      desc="uncorrected estimate of %s",
                      preferred="un",
                      ratio=ratio)
  }
  
  if(any(correction=="translate")) {
    # translation correction
    if(npts > 1) {
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=win, paired=TRUE)
      if(DEBUG) {
        elapsed <- proc.time() - started
        splat("pcf: computed edge weights after", codetime(elapsed))
      }
      kdenT <- sewpcf(d=dIJ, w=edgewt,
                         denargs=denargs,
                         lambda2area=lambda2area,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      if(DEBUG) {
        elapsed <- proc.time() - started
        splat("pcf: returned from sewpcf after", codetime(elapsed))
      }
      gT <- kdenT$g
      bw.used <- attr(kdenT, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenT, "bwvalues")
    } else gT <- undefined
    if(DEBUG) {
      elapsed <- proc.time() - started
      splat("Ready to bind.fv after", codetime(elapsed))
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(trans=gT),
                      denominator=samplesize,
                      labl="hat(%s)[Trans](r)",
                      desc="translation-corrected estimate of %s",
                      preferred="trans",
                      ratio=ratio)
  }

  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    if(npts > 1) {
      XI <- ppp(close$xi, close$yi, window=win, check=FALSE)
      edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
      if(DEBUG) {
        elapsed <- proc.time() - started
        splat("pcf: computed edge weights after", codetime(elapsed))
      }
      kdenR <- sewpcf(d=dIJ, w=edgewt,
                         denargs=denargs,
                         lambda2area=lambda2area,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      if(DEBUG) {
        elapsed <- proc.time() - started
        splat("pcf: returned from sewpcf after", codetime(elapsed))
      }
      gR <- kdenR$g
      bw.used <- attr(kdenR, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenR, "bwvalues")
    } else gR <- undefined
    if(DEBUG) {
      elapsed <- proc.time() - started
      splat("Ready to bind.fv after", codetime(elapsed))
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(iso=gR),
                      denominator=samplesize,
                      labl="hat(%s)[Ripley](r)",
                      desc="isotropic-corrected estimate of %s",
                      preferred="iso",
                      ratio=ratio)
  }
  
  if(DEBUG) {
    elapsed <- proc.time() - started
    splat("Finished bind.fv after", codetime(elapsed))
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }

  ## variance approximation
  ## Illian et al 2008 p 234 equation 4.3.42
  if(var.approx) {
    gr <- if(any(correction == "isotropic")) gR else gT
    # integral of squared kernel
    intk2 <- kernel.squint(kernel, bw.used)
    # isotropised set covariance of window
    gWbar <- as.function(rotmean(setcov(win), result="fv"))
    vest <- gr * intk2/(pi * r * gWbar(r) * lambda^2)
    out <- bind.ratfv(out,
                      quotient=data.frame(v=vest),
                      denominator=samplesize,
                      labl="v(r)", 
                      desc="approximate variance of %s",
                      ratio=ratio)
  }

  ## Finish off
  ## default is to display all corrections
  formula(out) <- . ~ r
  fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
  ##
  unitname(out) <- unitname(X)
  ## copy to other components
  if(ratio)
    out <- conform.ratfv(out)

  ## save information about computation
  attr(out, "bw") <- bw.used
  info <- append(info, list(bw.used=bw.used))
  if(adaptive) {
    attr(out, "bwvalues") <- bwvalues.used
    info <- append(info, list(bwvalues.used=bwvalues.used))
  }
  attr(out, "info") <- info
  if(DEBUG) {
    elapsed <- proc.time() - started
    splat("pcf: exiting, time taken", codetime(elapsed))
  }
  return(out)
}


# Smoothing Estimate of Weighted Pair Correlation
# d = vector of relevant distances
# w = vector of edge correction weights (in normal use)
# denargs = arguments to density.default (with defaults filled in)
# lambda2area = constant lambda^2 * areaW (in normal use)

sewpcf <- function(d, w, denargs, lambda2area,
                      divisor=c("r","d","a", "t"),
                      zerocor=c("none", "weighted", "convolution",
                                "reflection", "bdrykern", "JonesFoster"),
                      fast=TRUE, convert=TRUE, adaptive=FALSE, tau=0,
                      gref=NULL, Transform=NULL) {
  DEBUG <- isTRUE(getOption("debug.smoothpcf"))
  if(DEBUG) {
    started <- proc.time()
    splat("sewpcf")
  }

  ## trap outdated usage
  if(identical(as.character(divisor), c("r", "d"))) divisor <- "r" 
  divisor <- match.arg(divisor)
  zerocor <- match.arg(zerocor)
  w <- as.numeric(w)
  nw <- length(w)
  if(nw != length(d) && nw != 1)
    stop("Internal error: incorrect length of weights vector in sewpcf")
  if(divisor == "d") {
    w <- w/d
    if(!all(good <- is.finite(w))) {
      nbad <- sum(!good)
      warning(paste(nbad, "infinite, NA or NaN",
                    ngettext(nbad, "contribution was", "contributions were"),
                    "deleted from pcf estimate with divisor='d'.",
                    "Fraction deleted: ",
                    paste0(round(100 * nbad/length(w), 2), "%")),
              call.=FALSE)
      d <- d[good]
      w <- w[good]
    }
    nw <- length(w)
  }
  if(convert && divisor %in% c("a", "t") && is.numeric(denargs$bw)) {
    ## convert bandwidth value from 'distance' scale to 'area' scale
    switch(divisor,
           a = {
             denargs$bw <- with(denargs, pi * (from + to) * bw)
           },
           t = {
             midpoint <- with(denargs, (from + to)/2)
             midslope <- 2 * pi * midpoint * gref(midpoint)  # the integrand
             denargs$bw <- midslope * denargs$bw
           })
  }
  switch(divisor,
         r = {
           if(adaptive) {
             kden <- do.call(densityAdaptiveKernel,
                             append(list(quote(d),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           } else {
             kden <- do.call(densityBC,
                             append(list(quote(d),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           }
           r <- kden$x
           y <- kden$y
           if(tau != 0) # shrinkage
             y <- y + tau * 0.75/(kden$bw * sqrt(5))
           g <- y/(2 * pi * r * lambda2area)
         },
         d = {
           if(adaptive) {
             kden <- do.call(densityAdaptiveKernel,
                             append(list(quote(d),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           } else {
             kden <- do.call(densityBC,
                             append(list(quote(d),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           }
           r <- kden$x
           y <- kden$y 
           if(tau != 0) # shrinkage
             y <- y + tau * 0.75/(kden$bw * sqrt(5))
           g <- y/(2 * pi * lambda2area)
         },
         a = {
           areas <- pi * d^2
           ## interpret current density arguments
           r <- with(denargs, seq(from, to, length.out=n))
           ## adjust density arguments for computation
           if(convert) {
             denargs$n <- max(8192, 4 * denargs$n)
             denargs$to <- pi * denargs$to^2
             denargs$from <- pi * denargs$from^2
           }
           ## smooth on 'areas' scale
           if(adaptive) {
             kden <- do.call(densityAdaptiveKernel,
                             append(list(quote(areas),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           } else {
             kden <- do.call(densityBC,
                             append(list(quote(areas),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs)
                             )
           }
           ## linearly interpolate on 'areas' scale
           a.est <- kden$x
           g.est <- kden$y/lambda2area
           if(tau != 0) # shrinkage
             g.est <- g.est + tau * 0.75/(kden$bw * sqrt(5))/lambda2area
           a.query <- pi * r^2
           g <- approx(x = a.est,
                       y = g.est,
                       xout = a.query,
                       rule = c(2,1))$y
         },
         t = {
           ## transformation 
           tvalues <- Transform(d)
           ## interpret current density arguments
           r <- with(denargs, seq(from, to, length.out=n))
           ## adjust density arguments for computation
           if(convert) {
             denargs$n <- max(8192, 4 * denargs$n)
             denargs$to <- Transform(denargs$to)
             denargs$from <- Transform(denargs$from)
           }
           ## smooth on transformed scale
           if(adaptive) {
             kden <- do.call(densityAdaptiveKernel,
                             append(list(quote(tvalues),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs))
           } else {
             kden <- do.call(densityBC,
                             append(list(quote(tvalues),
                                         weights=quote(w),
                                         zerocor=zerocor,
                                         fast=fast),
                                    denargs))
           }
           ## linearly interpolate on transformed scale
           t.est <- kden$x
           y.est <- kden$y/lambda2area
           t.r <- Transform(r)
           y.r <- approx(x = t.est,
                         y = y.est,
                         xout = t.r,
                         rule = c(2,1))$y
           ## calculate g from y = g/gref                         
           g <- y.r * gref(r)
           if(tau != 0) # shrinkage
             g <- g + tau * 0.75/(denargs$bw * sqrt(5))/lambda2area
         })
  result <- data.frame(r=r,g=g)
  attr(result, "bw") <- kden$bw
  if(adaptive)
    attr(result, "bwvalues") <- kden$bwvalues
  if(DEBUG) {
    elapsed <- proc.time() - started
    splat("returning from sewpcf, time taken", codetime(elapsed))
  }
  return(result)
}

