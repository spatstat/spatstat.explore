#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.104 $	$Date: 2024/06/09 00:00:07 $
#
#	Kinhom()	compute estimate of K_inhom
#
#
#       Reference:
#            Non- and semiparametric estimation of interaction
#	     in inhomogeneous point patterns
#            A.Baddeley, J.Moller, R.Waagepetersen
#            Statistica Neerlandica 54 (2000) 329--350.
#
# -------- functions ----------------------------------------
#	Kinhom()	compute estimate of K
#                       using various edge corrections
#
#       Kwtsum()         internal routine for border correction
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
#       lambda          vector of intensity values for points of X
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
#			(denominator = sum of weights of points)
#
#       bord.modif:	K function estimated by border method
#			(denominator = area of eroded window)
#
# ------------------------------------------------------------------------

"Linhom" <- function(X, ..., correction) {
  if(missing(correction)) correction <- NULL
  K <- Kinhom(X, ..., correction=correction)
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, quote(L[inhom](r)), c("L", "inhom"),
                  names(K), new.labl=attr(K, "labl"))
  attr(L, "labl") <- attr(K, "labl")
  attr(L, "dangerous") <- attr(K, "dangerous")
  #
  return(L)  
}

"Kinhom"<-
  function (X, lambda=NULL, ..., r = NULL, breaks = NULL, 
            correction=c("border", "bord.modif", "isotropic", "translate"),
            renormalise=TRUE,
            normpower=1,
            update = TRUE,
            leaveoneout = TRUE,
            nlarge = 1000, 
            lambda2=NULL,
            reciplambda=NULL, reciplambda2=NULL,
	    diagonal=TRUE,
            sigma=NULL, varcov=NULL,
	    ratio=FALSE)
{
    verifyclass(X, "ppp")
    nlarge.given <- !missing(nlarge)
    rfixed <- !missing(r) || !missing(breaks)
    
    # determine basic parameters
    W <- X$window
    npts <- npoints(X)
    areaW <- area(W)
    diamW <- diameter(W)
    
    rmaxdefault <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    # match corrections
    correction.given <- !missing(correction) && !is.null(correction)
    if(is.null(correction))
      correction <- c("border", "bord.modif", "isotropic", "translate")

    correction <- pickoption("correction", correction,
                             c(none="none",
                               border="border",
                               "bord.modif"="bord.modif",
                               isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               good="good",
                               best="best"),
                             multi=TRUE)

#    best.wanted <- ("best" %in% correction)
    ## replace 'good' by the optimal choice for this size of dataset
    if("good" %in% correction)
      correction[correction == "good"] <- good.correction.K(X)
    ## retain only corrections that are implemented for the window
    correction <- implemented.for.K(correction, W$type, correction.given)

    ###########################################################
    # DETERMINE WEIGHTS AND VALIDATE
    #
    # The matrix 'lambda2' or 'reciplambda2' is sufficient information
    # unless we want the border correction.
    lambda2.given    <- !is.null(lambda2) || !is.null(reciplambda2)
    lambda2.suffices <- !any(correction %in% c("border", "bord.modif"))
    
    ## Arguments that are 'dangerous' for envelope, if fixed
    dangerous <- c("lambda", "reciplambda", "lambda2", "reciplambda2")
    danger <- TRUE
    
    # Use matrix of weights if it was provided and if it is sufficient
    if(lambda2.suffices && lambda2.given) {
      if(!is.null(reciplambda2)) {
        check.nmatrix(reciplambda2, npts, mname="reciplambda2")
        validate.weights(reciplambda2, recip=TRUE)
      } else {
        check.nmatrix(lambda2, npts, mname="lambda2")
        validate.weights(lambda2)
        reciplambda2 <- 1/lambda2
      }
      # renormalise
      if(renormalise && npts > 0) {
        check.1.real(normpower)
        stopifnot(normpower %in% 1:2)
	rlam2 <- reciplambda2
	if(!diagonal) diag(rlam2) <- 0
	renorm.factor <- (areaW^2/sum(rlam2))^(normpower/2)
      } 
    } else {
      ## Vector lambda or reciplambda is required
      a <- resolve.reciplambda(X, lambda=lambda, reciplambda=reciplambda,
                               ..., sigma=sigma, varcov=varcov,
                               leaveoneout=leaveoneout, update=update, check=TRUE)
      reciplambda <- a$reciplambda
      danger      <- a$danger
      dangerous   <- a$dangerous
      # renormalise
      if(renormalise && npts > 0) {
        check.1.real(normpower)
        stopifnot(normpower %in% 1:2)
        if(!diagonal && normpower == 2) {
	  renorm.factor <- (areaW^2)/(sum(reciplambda)^2 - sum(reciplambda^2))
	} else {
          renorm.factor <- (areaW/sum(reciplambda))^normpower
        }
      } 
    }

    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
        
  ###########################################
  # Efficient code for border correction and no correction
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

    can.do.fast <- breaks$even && !lambda2.given
    large.n    <- (npts >= nlarge)
#    demand.best <- correction.given && best.wanted
    large.n.trigger <- large.n && !correction.given
    fastcorrections <- c("border", "bord.modif", "none")
    fastdefault <- "border"
    correction.fast  <- all(correction %in% fastcorrections)
    will.do.fast <- can.do.fast && (correction.fast || large.n.trigger)
    asked.fast <- (correction.given && correction.fast) ||
                  (nlarge.given && large.n.trigger)
    if(!can.do.fast && asked.fast) {
      whynot <-
        if(!(breaks$even)) "r values not evenly spaced" else
        if(!missing(lambda)) "matrix lambda2 was given" else NULL
      warning(paste("cannot use efficient code", whynot, sep="; "))
    }
    if(will.do.fast) {
      ## Compute Kinhom using fast algorithm(s)
      ## determine correction(s)
      ok <- correction %in% fastcorrections
      correction <- if(any(ok)) correction[ok] else fastdefault
      bord <- any(correction %in% c("border", "bord.modif"))
      none <- any(correction =="none")
      if(!all(ok)) {
        ## some corrections were overridden; notify user
        corx <- c(if(bord) "border correction estimate" else NULL,
                  if(none) "uncorrected estimate" else NULL)
        corx <- paste(corx, collapse=" and ")
        message(paste("number of data points exceeds",
                      nlarge, "- computing", corx , "only"))
      }
      ## restrict r values to recommended range, unless specifically requested
      if(!rfixed) 
        r <- seq(from=0, to=alim[2], length.out=length(r))
      ## border method
      if(bord) {
        Kb <- Kborder.engine(X, max(r), length(r), correction,
                             weights=reciplambda, ratio=ratio)
        if(renormalise) {
          ynames <- setdiff(fvnames(Kb, "*"), "theo")
	  Kb <- adjust.ratfv(Kb, ynames, denfactor=1/renorm.factor)
        }
        Kb <- tweak.ratfv.entry(Kb, "border", new.labl="{hat(%s)[%s]^{bord}} (r)")
        Kb <- tweak.ratfv.entry(Kb, "bord.modif", new.labl="{hat(%s)[%s]^{bordm}} (r)")
      }
      ## uncorrected
      if(none) {
        Kn <- Knone.engine(X, max(r), length(r), weights=reciplambda,
	                   ratio=ratio)
        if(renormalise) 
	  Kn <- adjust.ratfv(Kn, "un", denfactor=1/renorm.factor)
        Kn <- tweak.ratfv.entry(Kn, "un", new.labl="{hat(%s)[%s]^{un}} (r)")
      }
      K <-
        if(bord && !none) Kb else
        if(!bord && none) Kn else
	if(!ratio) cbind.fv(Kb,  Kn[, c("r", "un")]) else 
	bind.ratfv(Kb,  Kn[, c("r", "un")], ratio=TRUE)
	
      ## tweak labels
      K <- rebadge.fv(K, quote(K[inhom](r)), c("K", "inhom"))
      if(danger)
        attr(K, "dangerous") <- dangerous
      return(K)
    }

  ###########################################
  # Fast code for rectangular window
  ###########################################

  if(can.do.fast && is.rectangle(W) && spatstat.options("use.Krect")) {
    K <-  Krect.engine(X, rmax, length(r), correction,
                        weights=reciplambda,
			ratio=ratio, fname=c("K", "inhom"))
    if(renormalise) {
      allfun <- setdiff(fvnames(K, "*"), "theo")
      K <- adjust.ratfv(K, allfun, denfactor=1/renorm.factor)
    }
    K <- rebadge.fv(K, quote(K[inhom](r)), c("K", "inhom"))
    attr(K, "alim") <- alim
    if(danger)
      attr(K, "dangerous") <- dangerous
    return(K)
  }
  
  ###########################################
  # Slower code
  ###########################################
        
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    denom <- if(renormalise) (areaW / renorm.factor) else areaW
    K <- ratfv(K, NULL, denom,
               argu="r",
	       ylab=quote(K[inhom](r)),
               valu="theo",
	       fmla=NULL,
	       alim=alim,
	       labl=c("r","{%s[%s]^{pois}}(r)"),
	       desc=desc,
               fname=c("K", "inhom"),
	       ratio=ratio)

    # identify all close pairs
    rmax <- max(r)
    what <- if(any(correction == "translate")) "all" else "ijd"
    close <- closepairs(X, rmax, what=what)
    dIJ <- close$d
    # compute weights for these pairs
    I <- close$i
    J <- close$j
#    wI <- reciplambda[I]
    wIJ <- 
      if(!lambda2.given)
        reciplambda[I] * reciplambda[J]
      else 
        reciplambda2[cbind(I,J)]
    # 

    # compute edge corrected estimates
    if(any(correction == "border" | correction == "bord.modif")) {
      # border method
      # Compute distances to boundary
      b <- bdist.points(X)
      bI <- b[I]
      # apply reduced sample algorithm
      RS <- Kwtsum(dIJ, bI, wIJ, b, w=reciplambda, breaks)
      if(any(correction == "border")) {
        Kb <- RS$ratio
        if(renormalise)
          Kb <- Kb * renorm.factor
        K <- bind.ratfv(K,
	                quotient = data.frame(border=Kb),
			denominator = denom,
	                labl = "{hat(%s)[%s]^{bord}}(r)",
                        desc = "border-corrected estimate of %s",
                        preferred = "border",
		        ratio=ratio)
      }
      if(any(correction == "bord.modif")) {
        Kbm <- RS$numerator/eroded.areas(W, r)
        if(renormalise)
          Kbm <- Kbm * renorm.factor
    	K <- bind.ratfv(K,
	                quotient = data.frame(bord.modif=Kbm),
			denominator = denom,
			labl = "{hat(%s)[%s]^{bordm}}(r)",
                        desc = "modified border-corrected estimate of %s",
                        preferred = "bord.modif",
			ratio=ratio)
      }
    }
    if(any(correction == "translate")) {
      # translation correction
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE)
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Ktrans <- cumsum(wh)/areaW
      if(renormalise)
        Ktrans <- Ktrans * renorm.factor
      rmax <- diamW/2
      Ktrans[r >= rmax] <- NA
      K <- bind.ratfv(K,
                      quotient = data.frame(trans=Ktrans),
		      denominator = denom,
		      labl ="{hat(%s)[%s]^{trans}}(r)",
                      desc = "translation-correction estimate of %s",
                      preferred = "trans",
		      ratio=ratio)
    }
    if(any(correction == "isotropic" | correction == "Ripley")) {
      # Ripley isotropic correction
      edgewt <- edge.Ripley(X[I], matrix(dIJ, ncol=1))
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Kiso <- cumsum(wh)/areaW
      if(renormalise)
        Kiso <- Kiso * renorm.factor
      rmax <- diamW/2
      Kiso[r >= rmax] <- NA
      K <- bind.ratfv(K,
                      quotient = data.frame(iso=Kiso),
		      denominator = denom,
		      labl = "{hat(%s)[%s]^{iso}}(r)",
                      desc = "Ripley isotropic correction estimate of %s",
                      preferred = "iso",
		      ratio=ratio)
    }

    # default is to display them all
    formula(K) <- . ~ r
    unitname(K) <- unitname(X)
    if(danger)
      attr(K, "dangerous") <- dangerous
    return(K)
}


Kwtsum <- function(dIJ, bI, wIJ, b, w, breaks, fatal=TRUE) {
  #
  # "internal" routine to compute border-correction estimates of Kinhom
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # wIJ:  product weight for selected I, J pairs
  #
  # b:    vector of ALL distances to window boundary
  # w:   weights for ALL points
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  stopifnot(length(bI) == length(wIJ))
  stopifnot(length(w) == length(b))

  if(!is.finite(sum(w, wIJ))) {
    if(fatal)
      stop("Weights in K-function were infinite or NA", call.=FALSE)
    #' set non-finite weights to zero
    if(any(bad <- !is.finite(w))) {
      warning(paste(sum(bad), "out of", length(bad),
                    paren(percentage(bad)), 
                    "of the boundary weights",
                    "in the K-function were NA or NaN or Inf",
                    "and were reset to zero"),
              call.=FALSE)
      w[bad] <- 0
    }
    if(any(bad <- !is.finite(wIJ))) {
      warning(paste(sum(bad), "out of", length(bad),
                    paren(percentage(bad)),
                    "of the weights for pairwise distances",
                    "in the K-function were NA or NaN or Inf",
                    "and were reset to zero"),
              call.=FALSE)
      wIJ[bad] <- 0
    }
  }
  
  bkval <- breaks$val
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  #
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], bkval, wIJ[uncen])
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], bkval, wIJ[uncen])
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, bkval, w)
  # total weight of censoring times beyond rightmost breakpoint
  uppercen <- sum(w[b > breaks$max])
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator   <- RS$numerator
  denominator <- RS$denominator
  ratio        <- RS$numerator/RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denominator) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  return(list(numerator=numerator, denominator=denominator, ratio=ratio))
}


