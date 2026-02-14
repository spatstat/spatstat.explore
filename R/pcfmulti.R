#
#   pcfmulti.R
#
#   $Revision: 1.16 $   $Date: 2026/02/14 05:50:13 $
#
#   multitype pair correlation functions
#

pcfcross <- function(X, i, j, ...) {
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))
  stopifnot(is.multitype(X))
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  ##
  result <- pcfmulti(X, I, J, ...,
                     Iname=Iname, Jname=Jname,
                     IJexclusive=(i != j))
  ##
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result,
               substitute(g[i,j](r),
                          list(i=iname,j=jname)),
               c("g", paste0("list", paren(paste(iname, jname, sep=",")))),
               new.yexp=substitute(g[list(i,j)](r),
                                   list(i=iname,j=jname)))
  return(result)
}

pcfdot <- 
function(X, i, ...) {
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))
  stopifnot(is.multitype(X))
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- "points"
	
  result <- pcfmulti(X, I, J, ...,
                     Iname=Iname, Jname=Jname,
                     IJexclusive=FALSE)

  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(g[i ~ dot](r), list(i=iname)),
               c("g", paste0(iname, "~symbol(\"\\267\")")),
               new.yexp=substitute(g[i ~ symbol("\267")](r),
                 list(i=iname)))
  return(result)
}


pcfmulti <- function(X, I, J, ...,
                     r=NULL, rmax=NULL,
                     adaptive=FALSE,
                     kernel="epanechnikov", bw=NULL, h=NULL,
                     bw.args=list(), stoyan=0.15, adjust=1,
                     correction=c("translate", "Ripley"),
                     divisor=c("a", "r", "d", "t"),
                     zerocor=c("convolution", "reflection", "bdrykern",
                               "JonesFoster", "weighted", "none",
                               "good", "best"),
                     nsmall = 300,
                     gref=NULL,
                     Iname="points satisfying condition I",
                     Jname="points satisfying condition J",
                     IJexclusive=FALSE,
                     ratio=FALSE)
{
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))

  win <- X$window
  areaW <- area(win)
  npts <- npoints(X)

  divisor.given <- !missing(divisor) && !is.null(divisor)
  zerocor.given <- !missing(zerocor) && !is.null(zerocor)
  correction.given <- !missing(correction) && !is.null(correction)

  if(!divisor.given || !zerocor.given) 
    warn.once("pcfmultiDefaults",
              paste("Default settings for pcfmulti, pcfdot, pcfcross",
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
  
  ## .......... indices I and J .............................
  
  I <- ppsubset(X, I, "I")
  J <- ppsubset(X, J, "J")
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")

  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  XI <- X[I]
  XJ <- X[J]

  lambdaJ <- nJ/areaW
  nIJ <- if(IJexclusive) 0 else sum(I & J)
  IJexclusive <- IJexclusive || (nIJ == 0)
  samplesize <- npairs <- nI * nJ - nIJ
  lambdaIJarea <- npairs/areaW
  

  ## ...........  kernel bandwidth and support .........................
  
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
                    bw <- stoyan/sqrt(5 * lambdaJ)
                    h <- bw * cker
                  },
                  bw.fiksel = ,
                  fiksel = {
                    ## Fiksel (1988)
                    bw <- 0.1/sqrt(lambdaJ)
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
  # handle argument r 

  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambdaJ)
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  ########## smoothing parameters for pcf ############################  
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
  hmax <- if(is.numeric(h)) h else (2*cker*stoyan/sqrt(lambdaJ))
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

  ## .......................................................
  ##  initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  fname <- c("g", "list(I,J)")
  yexp <- quote(g[list(I,J)](r))
  out <- ratfv(df=df, numer=NULL, denom=samplesize,
               argu="r",
               ylab=quote(g[I,J](r)),
               valu="theo", ,
               alim=alim,
               labl=c("r", makefvlabel(NULL, NULL, fname, "Pois")),
               desc=c("distance argument r", "theoretical Poisson %s"),
               fname=fname,
               yexp=yexp,
               ratio=ratio)
  
  #################################################
  
  ## compute pairwise distances
  
  ## identify close pairs of points
  if(npairs > 0) {
    what <- if(any(correction == "translate")) "all" else "ijd"
    if(IJexclusive) {
      close <- crosspairs(XI, XJ, rmax+hmax, what=what)
    } else {
      close <- crosspairs(XI, XJ, rmax+hmax, what=what,
                          iX=which(I), iY=which(J))
    }
    ## extract information for these pairs 
    dclose <- close$d
    icloseI  <- close$i  # (index in the sequence of XI)
  } else {
    undefined <- rep(NaN, length(r))
  }

  ###### compute #######

  bw.used <- bwvalues.used <- NULL

  if(any(correction=="none")) {
    ## uncorrected
    if(npairs > 0) {
      kdenN <- sewpcf(d=dclose, w=1, denargs, lambdaIJarea, divisor)
      gN <- kdenN$g
      bw.used <- attr(kdenN, "bw")
      if(adaptive) bwvalues.used <- attr(kdenN, "bwvalues")
    } else {
      gN <- undefined
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(un=gN),
                      denominator=samplesize,
                      labl=makefvlabel(NULL, "hat", fname, "un"),
                      desc="uncorrected estimate of %s",
                      preferred="un",
                      ratio=ratio)
  }
  if(any(correction=="translate")) {
    ## translation correction
    if(npairs > 0) {
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=win, paired=TRUE)
      kdenT <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)
      gT <- kdenT$g
      bw.used <- attr(kdenT, "bw")
      if(adaptive) bwvalues.used <- attr(kdenT, "bwvalues")
    } else {
      gT <- undefined
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(trans=gT),
                      denominator=samplesize,
                      labl=makefvlabel(NULL, "hat", fname, "Trans"),
                      desc="translation-corrected estimate of %s",
                      preferred="trans",
                      ratio=ratio)
  }
  if(any(correction=="isotropic")) {
    ## Ripley isotropic correction
    if(npairs > 0) {
      edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
      kdenR <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)
      gR <- kdenR$g
      bw.used <- attr(kdenR, "bw")
      if(adaptive) bwvalues.used <- attr(kdenR, "bwvalues")
    } else {
      gR <- undefined
    }
    out <- bind.ratfv(out,
                      quotient=data.frame(iso=gR),
                      denominator=samplesize,
                      labl=makefvlabel(NULL, "hat", fname, "Ripley"),
                      desc="isotropic-corrected estimate of %s",
                      preferred="iso",
                      ratio=ratio)
  }
  
  # which corrections have been computed?
  corrxns <- rev(setdiff(names(out), "r"))

  ## sanity check
  if(length(corrxns) == 0) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # default is to display them all
  formula(out) <- . ~ r
  fvnames(out, ".") <- corrxns

  # 
  unitname(out) <- unitname(X)
  return(out)
}

