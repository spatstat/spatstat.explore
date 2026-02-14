#
#   pcfmulti.R
#
#   $Revision: 1.17 $   $Date: 2026/02/14 09:04:30 $
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
                     tau = 0, 
                     Iname="points satisfying condition I",
                     Jname="points satisfying condition J",
                     IJexclusive=FALSE,
                     ratio=FALSE,
                     close=NULL)
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
  
  ## ..... distance values r ................................

  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambdaJ)
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  ## .... determine smoothing arguments ......................

  M <- resolve.pcf.bandwidth(X,
                             lambda=lambdaJ, 
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
    if(is.null(close)) {
      ## Find close pairs
      if(IJexclusive) {
        close <- crosspairs(XI, XJ, dmax, what=what)
      } else {
        close <- crosspairs(XI, XJ, dmax, what=what,
                            iX=which(I), iY=which(J))
      }
    } else {
      ## Use data provided
      needed <- if(what == "ijd") c("i", "j", "d") else
                 c("i", "j", "xi", "yi", "xj", "yj", "dx", "dy", "d")
      if(any(is.na(match(needed, names(close)))))
        stop(paste("Argument", sQuote("close"),
                   "should have components named",
                   commasep(sQuote(needed))),
             call.=FALSE)      
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
      kdenN <- sewpcf(d=dclose, w=1,
                      denargs=denargs,
                      lambda2area=lambdaIJarea,
                      divisor=divisor,
                      zerocor=zerocor,
                      adaptive=adaptive,
                      tau=tau,
                      gref=gref, Transform=Transform)
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
      kdenT <- sewpcf(d=dclose,
                      w=edgewt,
                      denargs=denargs,
                      lambda2area=lambdaIJarea,
                      divisor=divisor,
                      zerocor=zerocor,
                      adaptive=adaptive,
                      tau=tau,
                      gref=gref, Transform=Transform)
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
      kdenR <- sewpcf(d=dclose,
                      w=edgewt,
                      denargs=denargs,
                      lambda2area=lambdaIJarea,
                      divisor=divisor,
                      zerocor=zerocor,
                      adaptive=adaptive,
                      tau=tau,
                      gref=gref, Transform=Transform)
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

