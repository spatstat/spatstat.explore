#
#	Kest.R		Estimation of K function
#
#	$Revision: 5.139 $	$Date: 2022/06/30 07:49:47 $
#
#
# -------- functions ----------------------------------------
#	Kest()		compute estimate of K
#                       using various edge corrections
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
# -------- standard output ------------------------------
#      A data frame (class "fv") with columns named
#
#	r:		same as input
#
#	trans:		K function estimated by translation correction
#
#	iso:		K function estimated by Ripley isotropic correction
#
#	theo:		K function for Poisson ( = pi * r ^2 )
#
#	border:		K function estimated by border method
#			using standard formula (denominator = count of points)
#
#       bord.modif:	K function estimated by border method
#			using modified formula 
#			(denominator = area of eroded window
#
# ------------------------------------------------------------------------

"Lest" <- function(X, ..., correction) {
  if(missing(correction)) correction <- NULL
  K <- Kest(X, ..., correction=correction)
  L <- eval.fv(sqrt(K/pi), dotonly=FALSE)
  # handle variance estimates
  if(any(varcols <- colnames(K) %in% c("rip", "ls"))) {
    r <- with(L, .x)
    L[,varcols] <- as.data.frame(K)[,varcols]/(2 * pi * r)^2
    # fix 0/0
    n <- npoints(X)
    A <- area(Window(X))
    if(any(colnames(K) == "rip"))
      L[r == 0, "rip"] <- (2 * A/(n-1)^2)/(4 * pi)
    if(any(colnames(K) == "ls"))
      L[r == 0, "ls"]  <- (2 * A/(n * (n-1)))/(4 * pi)
  }
  # relabel the fv object
  L <- rebadge.fv(L, quote(L(r)), "L", names(K), new.labl=attr(K, "labl"))
  #
  return(L)  
}

"Kest"<-
function(X, ..., r=NULL, rmax=NULL, breaks=NULL, 
         correction=c("border", "isotropic", "Ripley", "translate"),
         nlarge=3000, domain=NULL, var.approx=FALSE,
         ratio=FALSE)
{
  verifyclass(X, "ppp")
  nlarge.given <- !missing(nlarge) && !is.null(nlarge)
  rfixed <- !is.null(r) || !is.null(breaks)
  npts <- npoints(X)
  npairs <-  npts * (npts - 1)
  W <- Window(X)
  areaW <- area(W)
  lambda <- npts/areaW
  lambda2 <- npairs/(areaW^2)
  lambda2area <- npairs/areaW
  samplesize <- npairs
  
  if(!is.null(domain)) {
    ## estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, W))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    ## use code in Kdot/Kmulti
    indom <- factor(inside.owin(X$x, X$y, domain), levels=c(FALSE,TRUE))
    Kd <- Kdot(X %mark% indom, i="TRUE",
               r=r, breaks=breaks, correction=correction,
               ratio=ratio, rmax=rmax,
               domainI=domain)
    # relabel and exit
    Kd <- rebadge.fv(Kd, quote(K(r)), "K")
    return(Kd)
  }

  rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambda)
  if(is.infinite(rmaxdefault)) rmaxdefault <- diameter(W)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  # choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("border", "isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             rigid="rigid",
                             periodic="periodic",
                             good="good",
                             best="best"),
                           multi=TRUE)
#  best.wanted <- ("best" %in% correction)
  # replace 'good' by the optimal choice for this size of dataset
  if("good" %in% correction)
    correction[correction == "good"] <- good.correction.K(X)
  # retain only corrections that are implemented for the window
  correction <- implemented.for.K(correction, W$type, correction.given)
  
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  ###########################################
  # Efficient code for border correction and no correction
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

  can.do.fast <- breaks$even
  large.n    <- (npts >= nlarge)
#  demand.best <- correction.given && best.wanted
  large.n.trigger <- large.n && !correction.given
  fastcorrections <- c("border", "bord.modif", "none")
  fastdefault     <- "border"
  correction.fast   <- all(correction %in% fastcorrections)
  will.do.fast <- can.do.fast && (correction.fast || large.n.trigger)
  asked <- correction.fast || (nlarge.given && large.n.trigger)
  if(asked && !can.do.fast)
    warning("r values not evenly spaced - cannot use efficient code")
  if(will.do.fast) {
    # determine correction(s)
    ok <- correction %in% fastcorrections
    correction <- if(any(ok)) correction[ok] else fastdefault
    bord <- any(correction %in% c("border", "bord.modif"))
    none <- any(correction =="none")
    if(!all(ok)) {
      # some corrections were overridden; notify user
      corx <- c(if(bord) "border correction estimate" else NULL,
                if(none) "uncorrected estimate" else NULL)
      corx <- paste(corx, collapse=" and ")
      message(paste("number of data points exceeds",
                    nlarge, "- computing", corx , "only"))
    }
    # restrict r values to recommended range, unless specifically requested
    if(!rfixed) 
      r <- seq(from=0, to=alim[2], length.out=length(r))
    if(bord)
      Kb <- Kborder.engine(X, max(r), length(r), correction, ratio=ratio)
    if(none)
      Kn <- Knone.engine(X, max(r), length(r), ratio=ratio)
    if(bord && none) {
      Kn <- Kn[ , names(Kn) != "theo"]
      yn <- fvnames(Kb, ".y")
      Kbn <- if(!ratio) bind.fv(Kb, Kn, preferred=yn) else
             bind.ratfv(Kb, Kn, preferred=yn)
      return(Kbn)
    }
    if(bord) return(Kb)
    if(none) return(Kn) 
  }

  unsupported.Krect <- c("rigid", "periodic")
  do.fast.rectangle <-
    can.do.fast &&
    is.rectangle(W) &&
    spatstat.options("use.Krect") &&
    !any(correction %in% unsupported.Krect)
  
  if(do.fast.rectangle) {
    ###########################################
    ## Fast code for rectangular window
    ###########################################
    K <-  Krect.engine(X, rmax, length(r), correction, ratio=ratio)
    attr(K, "alim") <- alim
  } else {
    ###########################################
    ## Slower code
    ###########################################

    ## this will be the output data frame
    Kdf <- data.frame(r=r, theo = pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    K <- ratfv(Kdf, NULL, npairs,
               "r", quote(K(r)),
               "theo", NULL, alim, c("r","%s[pois](r)"), desc, fname="K",
               ratio=ratio)
  
    ## Identify all close pairs up to distance 'rmax'
    rmax <- max(r)
    if(all(correction == "periodic")) {
      ## not needed in periodic case
      ## Assign null value to placate package checker
      close <- DIJ <- NULL
    } else {
      ## usual case
      ## Identify all close pairs
      needxy <- correction %in% c("translate", "isotropic")
      what <- if(any(needxy)) "all" else "ijd"
      close <- closepairs(X, rmax, what=what)
      DIJ <- close$d
    }
    
    ## precompute set covariance of window
    gW <- NULL
    if(any(correction %in% c("translate", "rigid", "isotropic")))
      gW <- setcov(W)
    
    if(any(correction == "none")) {
      ## uncorrected! For demonstration purposes only!
      wh <- whist(DIJ, breaks$val)  # no weights
      Kun <- cumsum(wh)/lambda2area
      ## uncorrected estimate of K
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(un=Kun),
                      denominator = npairs,
                      labl        = "hat(%s)[un](r)",
                      desc        = "uncorrected estimate of %s",
                      preferred   = "un",
                      ratio=ratio)
    }
  
    if(any(correction == "periodic")) {
      ## periodic correction
      ## Find close pairs of points in periodic distance
      closeP <- closepairs(X, rmax, what="ijd", periodic=TRUE)
      DIJP <- closeP$d
      ## Compute unweighted histogram
      wh <- whist(DIJP, breaks$val)
      Kper <- cumsum(wh)/lambda2area
      ## periodic correction estimate of K
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(per=Kper),
                      denominator = npairs,
                      labl        = "hat(%s)[per](r)",
                      desc        = "periodic-corrected estimate of %s",
                      preferred   = "per",
                      ratio=ratio)
    }
  
    if(any(correction == "border" | correction == "bord.modif")) {
      ## border method
      ## Compute distances to boundary
      b <- bdist.points(X)
      I <- close$i
      bI <- b[I]
      ## apply reduced sample algorithm
      RS <- Kount(DIJ, bI, b, breaks)
      if(any(correction == "bord.modif")) {
        ## modified border correction
        denom.area <- eroded.areas(W, r)
        Kbm <- RS$numerator/(lambda2 * denom.area)
        samplesizeKbm <- npairs * (denom.area/areaW)
        K <- bind.ratfv(K,
                        numerator   = NULL,
                        quotient    = data.frame(bord.modif=Kbm),
                        denominator = samplesizeKbm,
                        labl        = "hat(%s)[bordm](r)",
                        desc        = "modified border-corrected estimate of %s",
                        preferred   = "bord.modif",
                        ratio=ratio)
      }
      if(any(correction == "border")) {
        Kb <- RS$numerator/(lambda * RS$denom.count)
        samplesizeKb <- (npts-1) * RS$denom.count
        K <- bind.ratfv(K,
                        numerator   = NULL,
                        quotient    = data.frame(border=Kb), 
                        denominator = samplesizeKb,
                        labl        = "hat(%s)[bord](r)",
                        desc        = "border-corrected estimate of %s",
                        preferred   = "border",
                        ratio=ratio)
      }
    }

    if(any(correction == "translate")) {
      ## Ohser-Stoyan translation correction
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE,
                           gW = gW, give.rmax=TRUE)
      wh <- whist(DIJ, breaks$val, edgewt)
      Ktrans <- cumsum(wh)/lambda2area
      h <- attr(edgewt, "rmax")
      Ktrans[r >= h] <- NA
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(trans=Ktrans),
                      denominator = npairs,
                      labl        = "hat(%s)[trans](r)",
                      desc        = "translation-corrected estimate of %s",
                      preferred   = "trans",
                      ratio=ratio)
    }
    if(any(correction == "rigid")) {
      ## Ohser-Stoyan rigid motion correction
      CW <- rotmean(gW)
      edgewt <- areaW/as.function(CW)(DIJ)
      wh <- whist(DIJ, breaks$val, edgewt)
      Krigid <- cumsum(wh)/lambda2area
      h <- rmax.Rigid(X, gW) #sic: X not W
      Krigid[r >= h] <- NA
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(rigid=Krigid),
                      denominator = npairs,
                      labl        = "hat(%s)[rigid](r)",
                      desc        = "rigid motion-corrected estimate of %s",
                      preferred   = "rigid",
                      ratio=ratio)
    }
    if(any(correction == "isotropic")) {
      ## Ripley isotropic correction
      XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
      edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
      wh <- whist(DIJ, breaks$val, edgewt)
      Kiso <- cumsum(wh)/lambda2area
      h <- boundingradius(W)
      Kiso[r >= h] <- NA
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(iso=Kiso),
                      denominator = npairs,
                      labl        = "hat(%s)[iso](r)",
                      desc        = "Ripley isotropic correction estimate of %s",
                      preferred   = "iso",
                      ratio=ratio)
    }
  }

  #############################
  ##  VARIANCE APPROXIMATION
  #############################

  if(var.approx && !any(correction == "isotropic")) {
    warn.once("varapproxiso",
              "Ignored argument 'var.approx=TRUE'; the variance approximation",
              "is available only for the isotropic correction")
    var.approx <- FALSE
  }
  
  if(var.approx) {
    ## Compute variance approximations
    A <- areaW
    P <- perimeter(W)
    n <- npts
    ## Ripley asymptotic approximation
    rip <- 2 * ((A/(n-1))^2) * (pi * r^2/A + 0.96 * P * r^3/A^2
                                + 0.13 * (n/A) * P * r^5/A^2)
    ##
    vsamplesize <- (npts - 1)^2
    K <- bind.ratfv(K,
                    numerator   = NULL,
                    quotient    = data.frame(rip=rip),
                    denominator = vsamplesize,
                    labl        = "vR(r)", 
                    desc        = "Ripley approximation to var(%s) under CSR",
                    preferred   = "iso",
                    ratio=ratio)
    
    if(W$type == "rectangle") {
      ## Lotwick-Silverman
      a1r <- (0.21 * P * r^3 + 1.3 * r^4)/A^2
      a2r <- (0.24 * P * r^5 + 2.62 * r^6)/A^3
      ## contains correction to typo on p52 of Diggle 2003
      ## cf Lotwick & Silverman 1982 eq (5)
      br <- (pi * r^2/A) * (1 - pi * r^2/A) +
        (1.0716 * P * r^3 + 2.2375 * r^4)/A^2
      vls <- (A^2) * (2 * br - a1r + (n-2) * a2r)/(n*(n-1))
      ## add column 
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(ls=vls),
                      denominator = vsamplesize,
                      "vLS(r)",
                      "Lotwick-Silverman approx to var(%s) under CSR",
                      "iso",
                      ratio=ratio)
    }
  }
  
  ### FINISH OFF #####
  ## default plot will display all edge corrections
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  fvnames(K, ".") <- setdiff(nama, c("r", "rip", "ls"))
  ##
  unitname(K) <- unitname(X)
  # copy to other components
  if(ratio)
    K <- conform.ratfv(K)

  return(K)
}

################################################################  
#############  SUPPORTING ALGORITHMS ###########################
################################################################  

Kount <- function(dIJ, bI, b, breaks) {
  #
  # "internal" routine to compute border-correction estimate of K or Kij
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # b:    vector of ALL distances to window boundary
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], breaks$val)
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], breaks$val)
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, breaks$val)
  # count censoring times beyond rightmost breakpoint
  uppercen <- sum(b > max(breaks$val))
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator <- RS$numerator
  denom.count <- RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denom.count) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  
  return(list(numerator=numerator, denom.count=denom.count))
}

#### interface to C code for border method

Kborder.engine <- function(X, rmax, nr=100,
                           correction=c("border", "bord.modif"),
                           weights=NULL, ratio=FALSE) 
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  lambda <- npts/areaW
  npairs <- npts * (npts - 1)
  lambda2 <- npairs/(areaW^2)
  lambda2area <- npairs/areaW

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  Kfv <- ratfv(Kdf, NULL, npairs,
               "r", quote(K(r)),
               "theo",
               . ~ r,
               c(0,rmax), c("r","%s[pois](r)"), desc, fname="K",
               unitname=unitname(X),
               ratio=ratio)

  ####### start computing ############
  # sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  Xsort <- X[orderX]
  x <- Xsort$x
  y <- Xsort$y
  
  # boundary distances
  b <- bdist.points(Xsort)

  # call the C code
  if(is.null(weights)) {
    # determine whether the numerator can be stored as an integer
    bigint <- .Machine$integer.max
    if(npts < sqrt(bigint)) {
      # yes - use faster integer arithmetic
      res <- .C(SE_KborderI,
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                b=as.double(b),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.integer(integer(nr)),
                denom=as.integer(integer(nr)),
                PACKAGE="spatstat.explore")
    } else {
      # no - need double precision storage
      res <- .C(SE_KborderD,
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                b=as.double(b),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.double(numeric(nr)),
                denom=as.double(numeric(nr)),
                PACKAGE="spatstat.explore")
    }
    if("bord.modif" %in% correction) {
      denom.area <- eroded.areas(W, r)
      Kbm <- res$numer/(lambda2 * denom.area)
      samplesizeKbm <- npairs * (denom.area/areaW)
      Kfv <- bind.ratfv(Kfv,
                        numerator=NULL,
                        quotient=data.frame(bord.modif=Kbm),
                        denominator=samplesizeKbm,
                        "hat(%s)[bordm](r)",
                        "modified border-corrected estimate of %s",
                        "bord.modif",
                        ratio=ratio)
    }
    if("border" %in% correction) {
      Kb <- res$numer/(lambda * res$denom)
      samplesizeKb <- (npts - 1) * res$denom
      Kfv <- bind.ratfv(Kfv,
                        numerator=NULL,
                        quotient=data.frame(border=Kb),
                        denominator=samplesizeKb,
                        "hat(%s)[bord](r)",
                        "border-corrected estimate of %s",
                        "border",
                        ratio=ratio)
    }
  } else {
    ## weighted version
    if(is.numeric(weights)) {
      if(length(weights) != X$n)
        stop("length of weights argument does not match number of points in X")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(anyNA(weights))
        stop("domain of weights image does not contain all points of X")
    }
    weights.Xsort <- weights[orderX]
    res <- .C(SE_Kwborder,
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights.Xsort),
              b=as.double(b),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)),
              denom=as.double(numeric(nr)),
              PACKAGE="spatstat.explore")
    if("border" %in% correction) {
      numKb <- res$numer
      denKb <- res$denom
      Kfv <- bind.ratfv(Kfv,
                        numerator=data.frame(border=numKb),
                        denominator=data.frame(border=denKb),
                        "hat(%s)[bord](r)",
                        "border-corrected estimate of %s",
                        "border",
                        ratio=ratio)
    }
    if("bord.modif" %in% correction) {
      numKbm <- res$numer
      denKbm <- eroded.areas(W, r)
      Kfv <- bind.ratfv(Kfv,
                        numerator=data.frame(bord.modif=numKbm),
                        denominator=data.frame(bord.modif=denKbm),
                        "hat(%s)[bordm](r)",
                        "modified border-corrected estimate of %s",
                        "bord.modif",
                        ratio=ratio)
    }
  }
  return(Kfv)
}

Knone.engine <- function(X, rmax, nr=100,
                         weights=NULL, ratio=FALSE) 
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  ##  lambda <- npts/areaW
  npairs <- npts * (npts - 1)
  lambda2 <- npairs/(areaW^2)
  lambda2area <- npairs/areaW

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  Kfv <- ratfv(Kdf, NULL, npairs,
               "r", quote(K(r)),
               "theo",
               . ~ r,
               c(0,rmax), c("r","%s[pois](r)"), desc, fname="K",
               unitname=unitname(X),
               ratio=ratio)
  
  ####### start computing ############
  # sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  Xsort <- X[orderX]
  x <- Xsort$x
  y <- Xsort$y
  
  # call the C code
  if(is.null(weights)) {
    ## determine whether the numerator can be stored as an integer
    bigint <- .Machine$integer.max
    if(npts < sqrt(bigint)) {
      ## yes - use faster integer arithmetic
      res <- .C(SE_KnoneI,
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.integer(integer(nr)),
                PACKAGE="spatstat.explore")
    } else {
      ## no - need double precision storage
      res <- .C(SE_KnoneD,
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.double(numeric(nr)),
                PACKAGE="spatstat.explore")
    }
    Kun <- res$numer/lambda2area
    samplesizeKun <- npairs
  } else {
    ## weighted version
    if(is.numeric(weights)) {
      if(length(weights) != npts)
        stop("length of weights argument does not match number of points in X")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(anyNA(weights))
        stop("domain of weights image does not contain all points of X")
    }
    weights.Xsort <- weights[orderX]
    res <- .C(SE_Kwnone,
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights.Xsort),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)),
              PACKAGE="spatstat.explore")
    samplesizeKun <- totwt <- sum(weights)
    Kun <- res$numer/totwt
  }

  # tack on to fv object
  Kfv <- bind.ratfv(Kfv,
                    numerator=NULL,
                    quotient=data.frame(un=Kun),
                    denominator=samplesizeKun,
                    "hat(%s)[un](r)",
                    "uncorrected estimate of %s",
                    "un",
                    ratio=ratio)
  return(Kfv)
}

     

rmax.rule <- function(fun="K", W, lambda) {
  if(gotW <- !missing(W)) verifyclass(W, "owin")
  if(gotL <- !missing(lambda)) lambda <- as.numeric(lambda) # can be vector
  gotall <- gotW && gotL
  switch(fun,
         K = {
           ## Ripley's Rule
           ripley <- if(gotW) shortside(Frame(W))/4 else Inf
           ## Count at most 1000 neighbours per point
           rlarge <- if(gotL) sqrt(1000 /(pi * lambda)) else Inf
           rmax <- min(rlarge, ripley)
         },
         Kscaled = {
           ## rule of thumb for Kscaled
           rdiam  <- if(gotall) diameter(Frame(W))/2 * sqrt(lambda) else Inf
           rmax <- min(10, rdiam)
         },
         F = ,
         G = ,
         J = {
           # rule of thumb
           rdiam  <- if(gotW) diameter(Frame(W))/2 else Inf
           # Poisson process has F(rlarge) = 1 - 10^(-5)
           rlarge <- if(gotL) sqrt(log(1e5)/(pi * lambda)) else Inf
           rmax <- min(rlarge, rdiam)
         },
         stop(paste("Unrecognised function type", sQuote(fun)))
         )
  return(rmax)
}
           
    
implemented.for.K <- function(correction, windowtype, explicit) {
  if(any(b <- (correction == "best"))) {
    ## replace 'best' by the best available correction
    correction[b] <- switch(windowtype, mask="translate", "isotropic")
  }
  whinge <- NULL
  if(windowtype != "rectangle" && any(pe <- (correction == "periodic"))) {
    whinge <- "Periodic correction is not defined for non-rectangular windows"
    correction <- correction[!pe]
  }
  if(windowtype == "mask" && any(iso <- (correction == "iso"))) {
    whinge <- pasteN(whinge,
      "Isotropic correction is not implemented for binary mask windows",
      collapse=" and ")
    correction <- correction[!iso]
  }
  if(explicit && !is.null(whinge)) {
    if(length(correction)) {
      ## some desired corrections remain; warn about the deleted ones
      warning(whinge, call.=FALSE)
    } else {
      ## none of the desired corrections are supported
      stop(whinge, call.=FALSE)
    }
  }
  return(correction)
}

good.correction.K <- function(X) {
  nX <- npoints(X)
  W <- as.owin(X)
  avail <- c("none",
             if(nX < 1e5) "border" else NULL,
             if(nX < 3000)"translate" else NULL,
             if(nX < 1000 && !is.mask(W)) "isotropic" else NULL)
  chosen <- rev(avail)[1]
  return(chosen)
}

Krect.engine <- function(X, rmax, nr=100,
                         correction,
                         weights=NULL, ratio=FALSE, fname="K",
                         use.integers=TRUE) {
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  width <- sidelengths(W)[1]
  height <- sidelengths(W)[2]
  lambda <- npts/areaW
  npairs <- npts * (npts - 1)
  lambda2 <- npairs/(areaW^2)
  lambda2area <- npairs/areaW

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  if(weighted <- !is.null(weights)) {
    ## coerce weights to a vector
    if(is.numeric(weights)) {
      check.nvector(weights, npts, vname="weights")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(anyNA(weights))
        stop("domain of weights image does not contain all points of X")
    }
    totalweight <- sum(weights)
  }

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  denom <- if(weighted) areaW else (lambda2 * areaW)
  Kfv <- ratfv(Kdf, NULL, denom,
               "r", quote(K(r)),
               "theo",
               . ~ r,
               c(0,rmax),
               c("r", makefvlabel(NULL, NULL, fname, "pois")),
               desc, fname=fname,
               unitname=unitname(X),
               ratio=ratio)

  ####### prepare data ############

  if(!all(correction == "translate")) {
    ## Ensure rectangle has its bottom left corner at the origin
    if(W$xrange[1] != 0 || W$yrange[1] != 0) {
      X <- shift(X, origin="bottomleft")
      W <- as.owin(X)
    }
  }

  ## sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  x <- X$x[orderX]
  y <- X$y[orderX]
  if(weighted)
    wt <- weights[orderX]

  ## establish algorithm parameters
  doIso <- "isotropic" %in% correction 
  doTrans <- "translate" %in% correction
  doBord <- any(c("border", "bord.modif") %in% correction)
  doUnco <- "none" %in% correction
  trimedge <- spatstat.options("maxedgewt")

  ## allocate space for results
  ziso   <- numeric(if(doIso) nr else 1L)
  ztrans <- numeric(if(doTrans) nr else 1L)
  
  ## call the C code
  if(weighted) {
    ## weighted version
    zbnumer <- numeric(if(doBord) nr else 1L)
    zbdenom <- numeric(if(doBord) nr else 1L)
    zunco   <- numeric(if(doUnco) nr else 1L)
    res <- .C(SE_KrectWtd,
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(wt),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.double(zbnumer),
              bdenom=as.double(zbdenom),
              unco=as.double(zunco),
              PACKAGE="spatstat.explore")
  } else if(use.integers && npts < sqrt(.Machine$integer.max)) {
    ## unweighted
    ## numerator of border correction can be stored as an integer
    ## use faster integer arithmetic
    zbnumer <- integer(if(doBord) nr else 1L)
    zbdenom <- integer(if(doBord) nr else 1L)
    zunco   <- integer(if(doUnco) nr else 1L)
    res <- .C(SE_KrectInt,
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.integer(zbnumer),
              bdenom=as.integer(zbdenom),
              unco=as.integer(zunco),
              PACKAGE="spatstat.explore")
  } else {
    ## unweighted
    ## need double precision storage
    zbnumer <- numeric(if(doBord) nr else 1L)
    zbdenom <- numeric(if(doBord) nr else 1L)
    zunco   <- numeric(if(doUnco) nr else 1L)
    res <- .C(SE_KrectDbl,
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.double(zbnumer),
              bdenom=as.double(zbdenom),
              unco=as.double(zunco),
              PACKAGE="spatstat.explore")
  }

  ## Process corrections in reverse order of priority

  ## Uncorrected estimate
  if("none" %in% correction) {
    if(!weighted) {
      Kun <- res$unco/lambda2area
      samplesizeKun <- npairs
    } else {
      Kun <- res$unco/areaW
      samplesizeKun <- totalweight
    }
    Kfv <- bind.ratfv(Kfv,
                      numerator=NULL,
                      quotient = data.frame(un=Kun),
                      denominator=samplesizeKun,
                      makefvlabel(NULL, "hat", fname, "un"),
                      "uncorrected estimate of %s",
                      "un",
                      ratio=ratio)
  }
  
  ## Modified border correction
  if("bord.modif" %in% correction) {
    denom.area <- eroded.areas(W, r)
    if(!weighted) {
      Kbm <- res$bnumer/lambda2area
      samplesizeKbm <- npairs * (denom.area/areaW)
    } else {
      Kbm <- res$bnumer/denom.area
      samplesizeKbm <- denom.area
    }
    Kfv <- bind.ratfv(Kfv,
                      numerator=NULL,
                      quotient=data.frame(bord.modif=Kbm),
                      denominator=samplesizeKbm,
                      makefvlabel(NULL, "hat", fname, "bordm"),
                      "modified border-corrected estimate of %s",
                      "bord.modif",
                      ratio=ratio)
  }
  ## Border correction
  if("border" %in% correction) {
    if(!weighted) {
      Kb <- res$bnumer/(lambda * res$bdenom)
      samplesizeKb <- (npts - 1) * res$bdenom
    } else {
      Kb <- res$bnumer/res$bdenom
      samplesizeKb <- res$bdenom
    }
    Kfv <- bind.ratfv(Kfv,
                      numerator=NULL,
                      quotient=data.frame(border=Kb),
                      denominator=samplesizeKb,
                      makefvlabel(NULL, "hat", fname, "bord"),
                      "border-corrected estimate of %s",
                      "border",
                      ratio=ratio)
  }
  
  ## translation correction
  if("translate" %in% correction) {
    if(!weighted) {
      Ktrans <- res$trans/lambda2area
      samplesizeKtrans <- npairs
    } else {
      Ktrans <- res$trans/areaW
      samplesizeKtrans <- areaW
    }
    h <- diameter(as.rectangle(W))/2
    Ktrans[r >= h] <- NA
    Kfv <- bind.ratfv(Kfv,
                      numerator=NULL,
                      quotient=data.frame(trans=Ktrans),
                      denominator=samplesizeKtrans,
                      makefvlabel(NULL, "hat", fname, "trans"),
                      "translation-corrected estimate of %s",
                      "trans",
                      ratio=ratio)
  }
  ## isotropic correction
  if("isotropic" %in% correction) {
    if(!weighted) {
      Kiso <- res$iso/lambda2area
      samplesizeKiso <- npairs
    } else {
      Kiso <- res$iso/areaW
      samplesizeKiso <- areaW
    }
    h <- diameter(as.rectangle(W))/2
    Kiso[r >= h] <- NA
    Kfv <- bind.ratfv(Kfv,
                      numerator=NULL,
                      quotient=data.frame(iso=Kiso),
                      denominator=samplesizeKiso,
                      makefvlabel(NULL, "hat", fname, "iso"),
                      "isotropic-corrected estimate of %s",
                      "iso",
                      ratio=ratio)
  }
  ##
  return(Kfv)
}


  
