#'
#' pcfinhom.R
#'
#' Inhomogeneous pair correlation function of point pattern 
#'
#' Copyright (c) 2008-2025 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#' $Revision: 1.36 $ $Date: 2026/02/14 10:43:24 $

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
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))
  npts <- npoints(X)
  win <- Window(X)
  areaW <- area(win)
  
  kernel <- match.kernel(kernel)

  divisor.given <- !missing(divisor) && !is.null(divisor)
  zerocor.given <- !missing(zerocor) && !is.null(zerocor)
  correction.given <- !missing(correction) && !is.null(correction)

  if(!divisor.given || !zerocor.given) 
    warn.once("pcfinhomDefaults",
              paste("Default settings for pcfinhom",
                    "have changed in spatstat.explore 3.7-0.007"))

  if(divisor.given) {
    if(is.function(divisor)) divisor <- divisor(X)
    divisor <- match.arg(divisor)
  } else {
    divisor <- "a"
  }

  if(zerocor.given) {
    zerocor <- match.arg(zerocor)
    if(zerocor == "best") zerocor <- "JonesFoster"
    if(zerocor == "good") zerocor <- "convolution"
  } else {
    ## default depends on number of data points
    if(!missing(nsmall)) check.1.integer(nsmall)
    zerocor <- if(npts <= nsmall) "JonesFoster" else "convolution"
  }

  check.1.real(adjust)

  ## ..................................................
  ## ....... INTENSITY VALUES .........................
  ## ..................................................

  a <- resolve.reciplambda(X, lambda=lambda, reciplambda=reciplambda,
                           ..., sigma=sigma, adjust=adjust.sigma, varcov=varcov,
                           leaveoneout=leaveoneout, update=update, check=TRUE)
  reciplambda <- a$reciplambda
  lambda      <- a$lambda
  danger      <- a$danger
  dangerous   <- a$dangerous
  
  ## renormalise?
  if(renormalise && npts > 1) {
    if(missing(normpower)) {
      warn.once("pcfinhom.normpower",
                "Default value of normpower has changed in pcfinhom")
    } else {
      check.1.real(normpower)
      stopifnot(normpower %in% 1:2)
    }
    renorm.factor <- (areaW/sum(reciplambda))^normpower
  } 
  
  ## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  ## >>>>>>>>>>>>>> SPECIAL CASE <<<<<<<<<<<<<<<<<<<<<<<<<<
  ## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ## ....... handle argument 'domain' .......................  
  if(!is.null(domain)) {
    message("Sorry, argument 'domain' is currently not supported by pcfinhom")
    return(NAobject("fv")) ## CURRENTLY BROKEN
    ## estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, Window(X)))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    # trick pcfdot() into doing it
    indom <- inside.owin(X$x, X$y, domain)
    marx <- factor(indom, levels=c(FALSE,TRUE))
    g <- pcfdot.inhom(X %mark% marx,
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
                 "theo", NULL, c(0, max(g$r)), 
                 attr(g, "labl"), attr(g, "desc"), fname="g",
                 ratio=TRUE)
    }
    unitname(g) <- unitname(X)
    if(var.approx)
      warning("var.approx is not implemented when 'domain' is given")
    return(g)
  }

  ## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  ## >>>>>>>>>>>>>> NORMAL CASE <<<<<<<<<<<<<<<<<<<<<<<<<<<
  ## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  ## ...... get point pattern information .......
  verifyclass(X, "ppp")
  lambdaBar <- npts/areaW
  samplesize <- npairs <- npts * (npts - 1)
  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambdaBar)        

  ## ......... distance values .........................
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  ## ......... edge correction .........................
  if(correction.given) {
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
  } else {
    correction <- c("translate", "Ripley")
  }

  correction <- implemented.for.K(correction, win$type, correction.given)

  ## .... determine smoothing argyments ......................
  
  M <- resolve.pcf.bandwidth(X,
                             lambda=lambdaBar, ## NB
                             rmax=rmax, nr=length(r),
                             adaptive=adaptive, kernel=kernel,
                             bw=bw, h=h, bw.args=bw.args,
                             stoyan=stoyan, adjust=adjust,
                             correction=correction,
                             divisor=divisor,
                             zerocor=zerocor,
                             nsmall=nsmall,
                             gref=gref,
                             close=close)

  info    <- M$info
  denargs <- M$denargs

  Transform <- info$Transform
  dmax      <- info$dmax
  gref      <- info$gref

  
  #######################################################
  ## compute pairwise distances up to 'dmax'
  if(npairs > 0) {
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
    if(npairs > 0) {
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
    if(npairs > 0) {
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
    if(npairs > 0) {
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
  
  # which corrections have been computed?
  corrxns <- rev(setdiff(names(out), "r"))

  # sanity check
  if(length(corrxns) == 0) {
    warning("Nothing computed - no edge corrections chosen")
    return(NAobject("fv"))
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


