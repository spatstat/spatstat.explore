#
#   pcfmulti.inhom.R
#
#   $Revision: 1.23 $   $Date: 2026/02/14 11:06:54 $
#
#   inhomogeneous multitype pair correlation functions
#
#

pcfcross.inhom <- 
  function(X, i, j, lambdaI=NULL, lambdaJ=NULL, ...,
         r=NULL, breaks=NULL, rmax=NULL,
         kernel="epanechnikov", bw=NULL, adjust.bw=1, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         sigma=NULL, adjust.sigma=1, varcov=NULL)
{
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  g <- pcfmulti.inhom(X, I, J, lambdaI, lambdaJ, ...,
                      r=r,breaks=breaks, rmax=rmax,
                      kernel=kernel, bw=bw, adjust.bw=adjust.bw, stoyan=stoyan,
                      correction=correction,
                      sigma=sigma, adjust.sigma=adjust.sigma, varcov=varcov,
                      Iname=Iname, Jname=Jname)
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(g,
               substitute(g[inhom,i,j](r),
                          list(i=iname,j=jname)),
               c("g", paste0("list", paren(paste("inhom", i, j, sep=",")))),
               new.yexp=substitute(g[list(inhom,i,j)](r),
                                   list(i=iname,j=jname)))
  attr(result, "dangerous") <- attr(g, "dangerous")
  return(result)
}

pcfdot.inhom <- 
function(X, i, lambdaI=NULL, lambdadot=NULL, ...,
         r=NULL, breaks=NULL, rmax=NULL,
         kernel="epanechnikov", bw=NULL, adjust.bw=1, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         sigma=NULL, adjust.sigma=1, varcov=NULL)
{
  verifyclass(X, "ppp")
  if(is.NAobject(X)) return(NAobject("fv"))
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points")
	
  g <- pcfmulti.inhom(X, I, J, lambdaI, lambdadot, ...,
                      r=r,breaks=breaks, rmax=rmax,
                      kernel=kernel, bw=bw, adjust.bw=adjust.bw, stoyan=stoyan,
                      correction=correction,
                      sigma=sigma, adjust.sigma=adjust.sigma, varcov=varcov,
                      Iname=Iname, Jname=Jname)
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(g,
               substitute(g[inhom, i ~ dot](r), list(i=iname)),
               c("g", paste0("list(inhom,", iname, "~symbol(\"\\267\"))")),
               new.yexp=substitute(g[list(inhom, i ~ symbol("\267"))](r),
                 list(i=iname)))
  if(!is.null(dang <- attr(g, "dangerous"))) {
    dang[dang == "lambdaJ"] <- "lambdadot"
    dang[dang == "lambdaIJ"] <- "lambdaIdot"
    attr(result, "dangerous") <- dang
  }
  return(result)
}


pcfmulti.inhom <- function(X, I, J, lambdaI=NULL, lambdaJ=NULL, ...,
                           lambdaX=NULL,
                           r=NULL, breaks=NULL,  rmax=NULL,
                           adaptive=FALSE, 
                           kernel="epanechnikov",
                           bw=NULL, h=NULL, bw.args=list(),
                           stoyan=0.15, adjust.bw=1,
                           correction=c("translate", "Ripley"),
                           divisor=c("a", "r", "d", "t"),
                           zerocor=c("convolution", "reflection", "bdrykern",
                                     "JonesFoster", "weighted", "none",
                                     "good", "best"),
                           nsmall = 300,
                           sigma=NULL, adjust.sigma=1, varcov=NULL,
                           update=TRUE, leaveoneout=TRUE,
                           gref=NULL,
                           tau = 0,
                           Iname="points satisfying condition I",
                           Jname="points satisfying condition J",
                           IJexclusive=FALSE,
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

  check.1.real(adjust.bw)

  ##########  indices I and J  ########################
  
  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != npts || length(J) != npts)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  XI <- X[I]
  XJ <- X[J]

  nIJ <- if(IJexclusive) 0 else sum(I & J)
  npairs <- nI * nJ - nIJ
  IJexclusive <- IJexclusive || (nIJ == 0)

  lambdaJbar <- nJ/areaW
  
  ########## intensity values #########################
  a <- resolve.lambdacross(X=X, I=I, J=J,
                           lambdaI=lambdaI,
                           lambdaJ=lambdaJ,
                           lambdaX=lambdaX,
                           ...,
                           sigma=sigma, adjust=adjust.sigma, varcov=varcov,
                           leaveoneout=leaveoneout, update=update,
                           Iexplain=Iname, Jexplain=Jname)
  lambdaI <- a$lambdaI
  lambdaJ <- a$lambdaJ
  danger    <- a$danger
  dangerous <- a$dangerous

  ## .........distance values r ........................
  rmaxdefault <- rmax %orifnull% rmax.rule("K", win, lambdaJbar)        
  breaks <- handle.r.b.args(r, breaks, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  ## .............. edge correction .....................
  if(correction.given) {
    correction <- pickoption("correction", correction,
                             c(isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               good="translate",
                               best="best"),
                            multi=TRUE)
  } else {
    correction <- c("translate", "Ripley")
  }
  correction <- implemented.for.K(correction, win$type, correction.given)

  ## .... determine smoothing argyments ......................

  M <- resolve.pcf.bandwidth(X,
                             lambda=lambdaJbar, 
                             rmax=rmax, nr=length(r),
                             adaptive=adaptive, kernel=kernel,
                             bw=bw, h=h, bw.args=bw.args,
                             stoyan=stoyan, adjust=adjust.bw,
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

  
  #################################################
  
  ## initialise fv object
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  fname <- c("g", "list(inhom,I,J)")
  out <- fv(df, "r",
            quote(g[inhom,I,J](r)), "theo", ,
            alim,
            c("r", makefvlabel(NULL, NULL, fname, "pois")),            
            c("distance argument r", "theoretical Poisson %s"),
            fname=fname,
            yexp=quote(g[list(inhom,I,J)](r)))
  
  ## compute I-to-J distances up to 'dmax'

 if(npairs > 0) {
    needall <- any(correction %in% c("translate", "isotropic"))
    if(is.null(close)) {
      what <- if(needall) "all" else "ijd"
      if(IJexclusive) {
        close <- crosspairs(XI, XJ, dmax, what=what)
      } else {
        ## exclude pairs of identical points
        close <- crosspairs(XI, XJ, dmax, what=what,
                            iX=which(I), iY=which(J))
      }
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
    dclose <- close$d
    icloseI  <- close$i
    jcloseJ  <- close$j
  } else {
    undefined <- rep(NaN, length(r))
  }
  
  ## Form weight for each pair
  weight <- 1/(lambdaI[icloseI] * lambdaJ[jcloseJ])

  ###### compute #######

  bw.used <- bwvalues.used <- NULL  

  if(any(correction=="none")) {
    #' uncorrected
    if(npairs > 0) {
      kdenN <- sewpcf(dclose, edgewt * weight, denargs, areaW)
      gN <- kdenN$g
      bw.used <- attr(kdenN, "bw")
      if(adaptive) bwvalues.used <- attr(kdenN, "bwvalues")
    } else {
      gN <- undefined
    }
    out <- bind.fv(out,
                   data.frame(un=gN),
                   makefvlabel(NULL, "hat", fname, "un"),
                   "uncorrected estimate of %s",
                   "un")
  }
  if(any(correction=="translate")) {
    #' translation correction
    if(npairs > 0) {
      edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
      kdenT <- sewpcf(dclose, edgewt * weight, denargs, areaW)
      gT <- kdenT$g
      bw.used <- attr(kdenT, "bw")
      if(adaptive) bwvalues.used <- attr(kdenT, "bwvalues")
    } else {
      gT <- undefined
    }
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   makefvlabel(NULL, "hat", fname, "Trans"),
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    #' Ripley isotropic correction
    if(npairs > 0) {
      edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
      kdenR <- sewpcf(dclose, edgewt * weight, denargs, areaW)
      gR <- kdenR$g
      bw.used <- attr(kdenR, "bw")
      if(adaptive) bwvalues.used <- attr(kdenR, "bwvalues")
    } else {
      gR <- undefined
    }
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   makefvlabel(NULL, "hat", fname, "Ripley"),
                   "isotropic-corrected estimate of %s",
                   "iso")
  }
  
  # which corrections have been computed?
  corrxns <- rev(setdiff(names(out), "r"))

  # sanity check
  if(length(corrxns) == 0) {
    warning("Nothing computed - no edge corrections chosen")
    return(NAobject("fv"))
  }
  
  # default is to display them all
  formula(out) <- . ~ r
  fvnames(out, ".") <- corrxns

  #
  unitname(out) <- unitname(X)

  if(danger)
    attr(out, "dangerous") <- dangerous
  return(out)
}

