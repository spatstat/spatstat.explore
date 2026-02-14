#'
#' pcfinhom.R
#'
#' Inhomogeneous pair correlation function of point pattern 
#'
#' Copyright (c) 2008-2025 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#' $Revision: 1.34 $ $Date: 2026/02/14 01:51:37 $

pcfinhom <- function(X, lambda=NULL, ..., r=NULL, rmax=NULL, 
                     adaptive=FALSE,
                     kernel="epanechnikov", bw=NULL, h=NULL,
                     bw.args=list(),
                     stoyan=0.15,
                     adjust = 1,
                     correction=c("translate", "Ripley"),
                     divisor=c("a", "r", "d", "t"),
                     zerocor=c("convolution", "reflection", "bdrykern",
                               "JonesFoster", "weighted", "none",
                               "good", "best"),
                     nsmall = 300,
                     renormalise=TRUE,
                     normpower=2,
                     update=TRUE, leaveoneout=TRUE, reciplambda=NULL,
                     sigma=NULL, adjust.sigma=1, varcov=NULL,
                     gref=NULL,
                     tau = 0,
                     fast=TRUE,
                     var.approx=FALSE,
                     domain=NULL, ratio=FALSE,
                     close=NULL)
{
  if(is.NAobject(X)) return(NAobject("fv"))

  kernel <- match.kernel(kernel)
  if(is.function(divisor)) divisor <- divisor(X)
  if(divisor.given <- !missing(divisor)) {
    if(is.function(divisor)) divisor <- divisor(X)
    divisor <- match.arg(divisor)
  }
  if(zerocor.given <- !missing(zerocor)) {
    zerocor <- match.arg(zerocor)
    if(zerocor == "best") zerocor <- "JonesFoster"
    if(zerocor == "good") zerocor <- "convolution"
  }
  check.1.real(adjust)
  
  ## ...... get point pattern information .......
  verifyclass(X, "ppp")
  win <- Window(X)
  areaW <- area(win)
  npts <- npoints(X)
  lambdaBar <- npts/areaW
  samplesize <- npairs <- npts * (npts - 1)
  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambdaBar)        

  ## ....... handle argument 'domain' .......................  
  if(!is.null(domain)) {
    ## Apply different defaults in this case    
    if(!divisor.given) {
      divisor <- "d"
    } else if(!(divisor %in% c("r", "d"))) {
      stop(paste("Sorry, option divisor =", sQuote(divisor),
                 "is not yet available when 'domain' is given"),
           call.=FALSE)
    }
    if(!zerocor.given) {
      zerocor <- "none"
    } else if(zerocor != 'none') {
      stop(paste0("Sorry, option zerocor=", sQuote(zerocor),
                  "is not yet available when 'domain' is given"),
           call.=FALSE)
    }
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

  ## .............  Finally apply defaults for normal case .............
  
  if(!divisor.given || !zerocor.given) {

    warn.once("pcfinhomDefaults",
              "default settings have changed for pcfinhom",
              "in spatstat.explore >= 3.7-0.007")
    
    if(!divisor.given)
      divisor <- "a"

    if(!zerocor.given) {
      ## default depends on number of data points
      if(!missing(nsmall)) check.1.integer(nsmall)
      zerocor <- if(npts <= nsmall) "JonesFoster" else "convolution"
    }

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

  ## ....... intensity values .........................

  a <- resolve.reciplambda(X, lambda=lambda, reciplambda=reciplambda,
                           ..., sigma=sigma, adjust=adjust.sigma, varcov=varcov,
                           leaveoneout=leaveoneout, update=update, check=TRUE)
  reciplambda <- a$reciplambda
  lambda      <- a$lambda
  danger      <- a$danger
  dangerous   <- a$dangerous
  
  ## renormalise?
  if(renormalise && npts > 0) {
    if(missing(normpower)) {
      warn.once("pcfinhom.normpower",
                "Default value of normpower has changed in pcfinhom")
    } else {
      check.1.real(normpower)
      stopifnot(normpower %in% 1:2)
    }
    renorm.factor <- (areaW/sum(reciplambda))^normpower
  } 
  
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
                          lambda=lambda,
                          close=close)
             if(!("..." %in% bwformals)) 
               xtra <- xtra[intersect(names(xtra), bwformals)]
             bw.args <- resolve.defaults(bw.args, xtra)
             bw <- do.call(bw, append(list(quote(X)), bw.args))
             check.bandwidth(bw, "bandwidth value returned by function 'bw'")
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
                    bw <- stoyan/sqrt(5 * lambdaBar)
                    h <- bw * cker
                  },
                  bw.fiksel = ,
                  fiksel = {
                    ## Fiksel (1988)
                    bw <- 0.1/sqrt(lambdaBar)
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
                 stop("Internal error: pcfmodel() did not yield a function",
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
  hmax <- if(is.numeric(h) && length(h) == 1) h else
          (2*cker*stoyan/sqrt(lambdaBar))
  dmax <- rmax + hmax * if(kernel == "gaussian") 4 else 1
  if(is.numeric(denargs$bw)) {
    ## compute the bandwidth on the transformed scale now
    switch(divisor,
           r = {},
           d = {},
           a = {
             ## convert to bandwidth(s) for areas
             ## (using same rule as in 'sewpcf')
             bw.area <- with(denargs, pi * (from + to) * bw)
             dmax.area <- pi * rmax^2
             hmax.area <- cker * max(bw.area)
             dmax.area <- dmax.area + hmax.area
             dmax <- max(dmax, dmax.area)
           },
           t = {
             ## use transformation
             midpoint <- with(denargs, (from + to)/2)
             ## compute derivative of transformation at midpoint of interval
             midslope <- 2 * pi * midpoint * gref(midpoint)
             ## convert bandwidth to transformed scale
             bw.trans <- midslope * denargs$bw
             dmax.trans <-
               integrate(integrand, lower=0, upper=dmax, g=gref)$value
             if(!is.finite(dmax.trans))
               stop("Internal error: numerical quadrature failed",
                    call.=FALSE)
             hmax.trans <- cker * max(bw.trans)
             dmax.trans <- dmax.trans + hmax.trans
             dmax <- max(dmax, dmax.trans)
           })
  }
  info <- append(info, list(rmax=rmax, hmax=hmax, dmax=dmax))

  ## Precompute transform and inverse transform
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
    I <- close$i
    J <- close$j
    XI <- ppp(close$xi, close$yi, window=win, check=FALSE)
    wIJ <- reciplambda[I] * reciplambda[J]
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

  bw.used <- bwvalues.used <- NULL
  
  if(any(correction=="none")) {
    #' uncorrected
    if(npts > 1) {
      kdenN <- sewpcf(d=dIJ, w=wIJ,
                         denargs=denargs,
                         lambda2area=areaW,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      gN <- kdenN$g
      if(renormalise) gN <- gN * renorm.factor
      bw.used <- attr(kdenN, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenN, "bwvalues")
    } else gN <- undefined
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
      kdenT <- sewpcf(d=dIJ, w=edgewt * wIJ,
                         denargs=denargs,
                         lambda2area=areaW,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      gT <- kdenT$g
      if(renormalise) gT <- gT * renorm.factor
      bw.used <- attr(kdenT, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenT, "bwvalues")
    } else gT <- undefined
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
      kdenR <- sewpcf(d=dIJ, w=edgewt * wIJ,
                         denargs=denargs,
                         lambda2area=areaW,
                         divisor=divisor,
                         zerocor=zerocor,
                         fast=fast, adaptive=adaptive, tau=tau,
                         gref=gref, Transform=Transform)
      gR <- kdenR$g
      if(renormalise) gR <- gR * renorm.factor
      bw.used <- attr(kdenR, "bw")
      if(adaptive)
        bwvalues.used <- attr(kdenR, "bwvalues")
    } else gR <- undefined
    out <- bind.ratfv(out,
                      quotient=data.frame(iso=gR),
                      denominator=samplesize,
                      labl="hat(%s)[Ripley](r)",
                      desc="isotropic-corrected estimate of %s",
                      preferred="iso",
                      ratio=ratio)
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
    vest <- gr * intk2/(pi * r * gWbar(r) * lambdaBar^2)
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
  return(out)
}


