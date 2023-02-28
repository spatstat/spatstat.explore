#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 5.61 $	$Date: 2023/02/28 02:06:33 $
#
#
# -------- functions ----------------------------------------
#	Kcross()	cross-type K function K_{ij}
#                       between types i and j
#
#	Kdot()          K_{i\bullet}
#                       between type i and all points regardless of type
#
#       Kmulti()        (generic)
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#				including 'marks' vector
#	r		distance values at which to compute K	
#
# -------- standard output ------------------------------
#      A data frame with columns named
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

"Lcross" <- function(X, i, j, ..., from, to, correction) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- if(!missing(from)) from else levels(marks(X))[1]
  if(missing(j)) j <- if(!missing(to)) to else levels(marks(X))[2]
  if(missing(correction)) correction <- NULL
  K <- Kcross(X, i, j, ..., correction=correction)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  L <- rebadge.fv(L,
                  substitute(L[i,j](r),
                             list(i=iname,j=jname)),
                  c("L", paste0("list(", iname, ",", jname, ")")),
                  new.yexp=substitute(L[list(i,j)](r),
                                      list(i=iname,j=jname)))
  attr(L, "labl") <- attr(K, "labl")
  return(L)  
}

"Ldot" <- function(X, i, ..., from, correction) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- if(!missing(from)) from else levels(marks(X))[1]
  if(missing(correction)) correction <- NULL
  K <- Kdot(X, i, ..., correction=correction)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  iname <- make.parseable(paste(i))
  L <- rebadge.fv(L,
                  substitute(L[i ~ dot](r), list(i=iname)),
                  c("L", paste(iname, "~ symbol(\"\\267\")")), 
                  new.yexp=substitute(L[i ~ symbol("\267")](r), list(i=iname)))
  attr(L, "labl") <- attr(K, "labl")
  return(L)  
}

"Kcross" <- 
function(X, i, j, r=NULL, breaks=NULL,
         correction =c("border", "isotropic", "Ripley", "translate") , ...,
         ratio=FALSE, from, to)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- if(!missing(from)) from else levels(marx)[1]
  if(missing(j))
    j <- if(!missing(to)) to else levels(marx)[2]
  I <- (marx == i)
  if(!any(I))
    stop(paste("No points have mark i =", i))

  if(i == j) {
    ## use Kest
    XI <- X[I]
    dont.complain.about(XI)
    result <- do.call(Kest,
                      resolve.defaults(list(X=quote(XI),
                                            r=quote(r),
                                            breaks=quote(breaks),
                                            correction=correction, ratio=ratio),
                                       list(rmax=NULL), ## forbidden 
                                       list(...)))
  } else {
    J <- (marx == j)
    if(!any(J))
      stop(paste("No points have mark j =", j))
    result <- Kmulti(X, I, J,
                     r=r, breaks=breaks,
                     correction=correction, ratio=ratio, ...)
  }
  result <- rebadge.as.crossfun(result, "K", NULL, i, j)
  return(result)
}

"Kdot" <- 
function(X, i, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...,
         ratio=FALSE, from)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- if(!missing(from)) from else levels(marx)[1]
        
  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
	
  if(!any(I)) stop(paste("No points have mark i =", i))
	
  result <- Kmulti(X, I, J,
                   r=r, breaks=breaks, correction=correction, ..., ratio=ratio)
  result <- rebadge.as.dotfun(result, "K", NULL, i)
  return(result)
}


"Kmulti"<-
function(X, I, J, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...,
         rmax=NULL, ratio=FALSE)
{
  verifyclass(X, "ppp")

  npts <- npoints(X)
  W <- Window(X)
  areaW <- area(W)

  dotargs <- list(...)
  domainI <- resolve.1.default("domainI", dotargs) %orifnull% W
  domainJ <- resolve.1.default("domainJ", dotargs) %orifnull% W
  areaI <- area(domainI)
  areaJ <- area(domainJ)

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
                             periodic="periodic",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  I <- ppsubset(X, I, "I")
  J <- ppsubset(X, J, "J")
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
	
  if(!any(I)) stop("no points belong to subset I")
  if(!any(J)) stop("no points belong to subset J")
		
  nI <- sum(I)
  nJ <- sum(J)
  lambdaI <- nI/areaI
  lambdaJ <- nJ/areaJ
  npairs <- nI * nJ

  # r values 
  rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambdaJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  # It will be given more columns later
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- ratfv(K, NULL, npairs,
             "r", quote(K[IJ](r)), 
             "theo", , alim, c("r","{%s[%s]^{pois}}(r)"),
             desc, fname=c("K", "list(I,J)"),
             yexp=quote(K[list(I,J)](r)),
             ratio=ratio)
  
  ## Extract relevant points
  XI <- X[I]
  XJ <- X[J]
  ## Map XI and XJ to original serial numbers in X
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]

  ## Find close pairs of points
  if(any(correction != "periodic")) {
    ## Find close pairs of points in Euclidean distance
    close <- crosspairs(XI, XJ, max(r), what="ijd",
                        iX=imap, iY=jmap)
    ## extract information for these pairs (relative to orderings of XI, XJ)
    dcloseIJ <- close$d
    icloseI  <- close$i
    jcloseJ  <- close$j
  }

  ## ...........................................................
  ## Compute estimates by each of the selected edge corrections.
  ## ...........................................................

  if(any(correction == "none")) {
    # uncorrected! 
    wh <- whist(dcloseIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/(lambdaI * lambdaJ * areaI)
    K <- bind.ratfv(K,
                    numerator   = NULL,
                    quotient    = data.frame(un=Kun),
                    denominator = npairs,
                    labl        = "{hat(%s)[%s]^{un}}(r)",
                    desc        = "uncorrected estimate of %s",
                    preferred   = "un",
                    ratio=ratio)
  }
  
  if(any(correction == "periodic")) {
    ## Periodic (toroidal) correction
    ## Compute periodic distances
    closeP <- crosspairs(XI, XJ, max(r), what="ijd",
                         periodic=TRUE,
                         iX=imap, iY=jmap)
    ## evaluate estimate
    wh <- whist(closeP$d, breaks$val)  # no weights
    Kper <- cumsum(wh)/(lambdaI * lambdaJ * areaI)
    K <- bind.ratfv(K,
                    numerator   = NULL,
                    quotient    = data.frame(per=Kper),
                    denominator = npairs,
                    labl        = "{hat(%s)[%s]^{per}}(r)",
                    desc        = "periodic-corrected estimate of %s",
                    preferred   = "per",
                    ratio=ratio)
  }
  
  if(any(correction == "border" | correction == "bord.modif")) {
    # border method
    # distance to boundary from each point of type I
    bI <- bdist.points(XI)
    # distance to boundary from first element of each (i, j) pair
    bcloseI <- bI[icloseI]
    # apply reduced sample algorithm
    RS <- Kount(dcloseIJ, bcloseI, bI, breaks)
    if(any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      Kbm <- RS$numerator/(denom.area * npairs)
      samplesizeKbm <- npairs * (denom.area/areaW)
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(bord.modif=Kbm),
                      denominator = npairs,
                      labl        = "{hat(%s)[%s]^{bordm}}(r)",
                      desc        = "modified border-corrected estimate of %s",
                      preferred   = "bord.modif",
                      ratio=ratio)
    }
    if(any(correction == "border")) {
      Kb <- RS$numerator/(lambdaJ * RS$denom.count)
      K <- bind.ratfv(K,
                      numerator   = NULL,
                      quotient    = data.frame(border=Kb),
                      denominator = RS$denom.count * nJ,
                      labl        = "{hat(%s)[%s]^{bord}}(r)",
                      desc        = "border-corrected estimate of %s",
                      preferred   = "border",
                      ratio=ratio)
    }
  }
  if(any(correction == "translate")) {
    # translation correction
    edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/(lambdaI * lambdaJ * areaI)
    rmax <- diameter(W)/2
    Ktrans[r >= rmax] <- NA
    K <- bind.ratfv(K,
                    numerator   = NULL,
                    quotient    = data.frame(trans=Ktrans), 
                    denominator = npairs,
                    labl        = "{hat(%s)[%s]^{trans}}(r)",             
                    desc        = "translation-corrected estimate of %s", 
                    preferred   = "trans",                                
                    ratio=ratio)
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dcloseIJ, ncol=1))
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/(lambdaI * lambdaJ * areaI)
    rmax <- diameter(W)/2
    Kiso[r >= rmax] <- NA
    K <- bind.ratfv(K,
                    numerator   = NULL,
                    quotient    = data.frame(iso=Kiso),
                    denominator = npairs,
                    labl        = "{hat(%s)[%s]^{iso}}(r)",             
                    desc        = "Ripley isotropic corrected estimate of %s", 
                    preferred   = "iso",                                
                    ratio=ratio)
  }
  # default is to display them all
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  
  if(ratio) K <- conform.ratfv(K)
  
  return(K)
}


