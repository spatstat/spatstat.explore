#
#	Kmulti.inhom.S		
#
#	$Revision: 1.56 $	$Date: 2023/02/28 02:06:17 $
#
#
# ------------------------------------------------------------------------

Lcross.inhom <- function(X, i, j, ..., correction) {
  if(!is.multitype(X, dfok=FALSE))
	stop("Point pattern must be multitype")
  if(missing(i)) i <- levels(marks(X))[1]
  if(missing(j)) j <- levels(marks(X))[2]
  if(missing(correction)) correction <- NULL
  K <- Kcross.inhom(X, i, j, ..., correction=correction)
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  # relabel the fv object
  L <- rebadge.fv(L,
                  substitute(L[inhom,i,j](r),
                             list(i=iname,j=jname)),
                  c("L", paste0("list", paren(paste("inhom", i, j, sep=",")))),
                  new.yexp=substitute(L[list(inhom,i,j)](r),
                                      list(i=iname,j=jname)))
  attr(L, "labl") <- attr(K, "labl")
  attr(L, "dangerous") <- attr(K, "dangerous")
  return(L)  
}

Ldot.inhom <- function(X, i, ..., correction) {
  if(!is.multitype(X, dfok=FALSE))
	stop("Point pattern must be multitype")
  if(missing(i)) i <- levels(marks(X))[1]
  if(missing(correction)) correction <- NULL
  K <- Kdot.inhom(X, i, ..., correction=correction)
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  # relabel the fv object
  iname <- make.parseable(paste(i))
  L <- rebadge.fv(L,
                  substitute(L[inhom, i ~ dot](r), list(i=iname)),
                  c("L", paste0("list(inhom,", iname, "~symbol(\"\\267\"))")),
                  new.yexp=substitute(L[list(inhom, i ~ symbol("\267"))](r),
                    list(i=iname)))
  attr(L, "labl") <- attr(K, "labl")
  attr(L, "dangerous") <- attr(K, "dangerous")
  return(L)  
}

"Kcross.inhom" <- 
function(X, i, j, lambdaI=NULL, lambdaJ=NULL, ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate"),
         sigma=NULL, varcov=NULL,
         lambdaIJ=NULL,
         lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE))
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL
  miss.update <- missing(update)
  miss.leave <- missing(leaveoneout)
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  K <- Kmulti.inhom(X, I, J, lambdaI, lambdaJ, ...,
                    r=r,breaks=breaks,correction=correction,
                    sigma=sigma, varcov=varcov,
                    lambdaIJ=lambdaIJ, Iname=Iname, Jname=Jname,
                    lambdaX=lambdaX, update=update, leaveoneout=leaveoneout,
                    miss.update=miss.update, miss.leave=miss.leave)
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(K,
               substitute(K[inhom,i,j](r),
                          list(i=iname,j=jname)),
               c("K", paste0("list", paren(paste("inhom", i, j, sep=",")))),
               new.yexp=substitute(K[list(inhom,i,j)](r),
                                   list(i=iname,j=jname)))
  attr(result, "dangerous") <- attr(K, "dangerous")
  return(result)
}

"Kdot.inhom" <- 
function(X, i, lambdaI=NULL, lambdadot=NULL, ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate"),
         sigma=NULL, varcov=NULL, 
         lambdaIdot=NULL,
         lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE))
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL
  miss.update <- missing(update)
  miss.leave <- missing(leaveoneout)

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points")
	
  K <- Kmulti.inhom(X, I, J, lambdaI, lambdadot, ...,
                    r=r,breaks=breaks,correction=correction,
                    sigma=sigma, varcov=varcov,
                    lambdaIJ=lambdaIdot,
                    Iname=Iname, Jname=Jname,
                    lambdaX=lambdaX, update=update, leaveoneout=leaveoneout,
                    miss.update=miss.update, miss.leave=miss.leave)
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(K,
               substitute(K[inhom, i ~ dot](r), list(i=iname)),
               c("K", paste0("list(inhom,", iname, "~symbol(\"\\267\"))")),
               new.yexp=substitute(K[list(inhom, i ~ symbol("\267"))](r),
                                   list(i=iname)))
  if(!is.null(dang <- attr(K, "dangerous"))) {
    dang[dang == "lambdaJ"] <- "lambdadot"
    dang[dang == "lambdaIJ"] <- "lambdaIdot"
    attr(result, "dangerous") <- dang
  }
  return(result)
}


"Kmulti.inhom"<-
function(X, I, J, lambdaI=NULL, lambdaJ=NULL, 
         ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate"),
         lambdaIJ=NULL,
         sigma=NULL, varcov=NULL,
         lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
{
  verifyclass(X, "ppp")

  dflt <- list(Iname="points satisfying condition I",
               Jname="points satisfying condition J",
               miss.update=missing(update),
               miss.leave=missing(leaveoneout))

  extrargs <- resolve.defaults(list(...), dflt)
  if(length(extrargs) > length(dflt))
    warning("Additional arguments unrecognised")
  Iname <- extrargs$Iname
  Jname <- extrargs$Jname
  miss.update <- extrargs$miss.update
  miss.leave <- extrargs$miss.leave
        
  npts <- npoints(X)
  W <- as.owin(X)
  areaW <- area(W)

  # validate edge correction
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
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  # validate I, J
  I <- ppsubset(X, I, "I")
  J <- ppsubset(X, J, "J")
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
  XI <- X[I]
  XJ <- X[J]
  
  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  ## r values 
  rmaxdefault <- rmax.rule("K", W, nJ/areaW)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  ## lambda values
  a <- resolve.lambdacross(X=X, I=I, J=J,
                           lambdaI=lambdaI,
                           lambdaJ=lambdaJ,
                           lambdaX=lambdaX,
                           lambdaIJ=lambdaIJ,
                           ...,
                           sigma=sigma, varcov=varcov,
                           leaveoneout=leaveoneout, update=update,
                           Iexplain=Iname, Jexplain=Jname)
  lambdaI <- a$lambdaI
  lambdaJ <- a$lambdaJ
  lambdaIJ <- a$lambdaIJ
  danger    <- a$danger
  dangerous <- a$dangerous

  ## Weight for each pair
  if(!is.null(lambdaIJ)) {
    danger <- TRUE
    dangerous <- union(dangerous, "lambdaIJ")
    if(!is.matrix(lambdaIJ))
      stop("lambdaIJ should be a matrix")
    if(nrow(lambdaIJ) != nI)
      stop(paste("nrow(lambdaIJ) should equal the number of", Iname))
    if(ncol(lambdaIJ) != nJ)
      stop(paste("ncol(lambdaIJ) should equal the number of", Jname))
  }

  # Recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  # It will be given more columns later
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  fname <- c("K", "list(inhom,I,J)")
  K <- fv(K, "r", quote(K[inhom, I, J](r)),
          "theo", , alim,
          c("r", makefvlabel(NULL, NULL, fname, "pois")),
          desc,
          fname=fname,
          yexp=quote(K[list(inhom,I,J)](r)))

# identify close pairs of points
  close <- crosspairs(XI, XJ, max(r), what="ijd")
# map (i,j) to original serial numbers in X
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
# eliminate any identical pairs
  if(any(I & J)) {
    ok <- (iX != jX)
    if(!all(ok)) {
      close$i  <- close$i[ok]
      close$j  <- close$j[ok]
      close$d  <- close$d[ok]
    }
  }
# extract information for these pairs (relative to orderings of XI, XJ)
  dclose <- close$d
  icloseI  <- close$i
  jcloseJ  <- close$j
        
# Form weight for each pair
  if(is.null(lambdaIJ))
    weight <- 1/(lambdaI[icloseI] * lambdaJ[jcloseJ])
  else 
    weight <- 1/lambdaIJ[cbind(icloseI, jcloseJ)]

# Compute estimates by each of the selected edge corrections.

  if(any(correction == "none")) {
    ## uncorrected
    wh <- whist(dclose, breaks$val, weight)
    Kun <- cumsum(wh)/areaW
    rmax <- diameter(W)/2
    Kun[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(un=Kun),
                 makefvlabel(NULL, "hat", fname, "un"),
                 "uncorrected estimate of %s",
                 "un")
  }
  
  if(any(correction == "border" | correction == "bord.modif")) {
    # border method
    # Compute distances to boundary
    b <- bdist.points(XI)
    bI <- b[icloseI]
    # apply reduced sample algorithm
    RS <- Kwtsum(dclose, bI, weight, b, 1/lambdaI, breaks)
    if(any(correction == "border")) {
      Kb <- RS$ratio
      K <- bind.fv(K, data.frame(border=Kb),
                   makefvlabel(NULL, "hat", fname, "bord"),
                   "border-corrected estimate of %s",
                   "border")
    }
    if(any(correction == "bord.modif")) {
      Kbm <- RS$numerator/eroded.areas(W, r)
      K <- bind.fv(K, data.frame(bord.modif=Kbm),
                   makefvlabel(NULL, "hat", fname, "bordm"),
                   "modified border-corrected estimate of %s",
                   "bord.modif")
    }
  }
  if(any(correction == "translate")) {
    ## translation correction
    edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
    allweight <- edgewt * weight
    wh <- whist(dclose, breaks$val, allweight)
    Ktrans <- cumsum(wh)/areaW
    rmax <- diameter(W)/2
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans),
                 makefvlabel(NULL, "hat", fname, "trans"),
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    ## Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
    allweight <- edgewt * weight
    wh <- whist(dclose, breaks$val, allweight)
    Kiso <- cumsum(wh)/areaW
    rmax <- diameter(W)/2
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), 
                 makefvlabel(NULL, "hat", fname, "iso"),
                 "Ripley isotropic correction estimate of %s",
                 "iso")
  }
  ## default is to display them all
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  if(danger)
    attr(K, "dangerous") <- dangerous
  return(K)
}
