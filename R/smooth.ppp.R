#
#  smooth.ppp.R
#
#  Smooth the marks of a point pattern
# 
#  $Revision: 1.90 $  $Date: 2024/06/09 00:02:42 $
#

Smooth <- function(X, ...) {
  UseMethod("Smooth")
}

Smooth.solist <- function(X, ...) {
  solapply(X, Smooth, ...)
}

Smooth.ppp <- function(X, sigma=NULL, ...,
                       weights=rep(1, npoints(X)), at="pixels",
                       leaveoneout=TRUE,
                       adjust=1, varcov=NULL, 
                       edge=TRUE, diggle=FALSE, 
                       kernel="gaussian",
                       scalekernel=is.character(kernel),
                       se=FALSE,
                       loctype=c("random", "fixed"),
                       wtype=c("multiplicity", "importance"),
                       geometric=FALSE) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE, na.action="fatal"))
    stop("X should be a marked point pattern", call.=FALSE)
  nX <- npoints(X)
  X <- coerce.marks.numeric(X)
  marx <- marks(X)
  if(!all(is.finite(as.matrix(marx))))
    stop("Some mark values are Inf, NaN or NA", call.=FALSE)
  univariate <- is.null(dim(marx))

  ## options
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  loctype <- match.arg(loctype)
  wtype <- match.arg(wtype)

  ## trivial case
  if(nX == 0) {
    cn <- colnames(marks(X))
    nc <- length(cn)
    switch(at,
           points = {
             Estimate <- if(univariate) numeric(0) else
                       matrix(, 0, nc, dimnames=list(NULL, cn))
             result <- if(!se) Estimate else
                       cbind(estimate=Estimate, SE=Estimate)
           },
           pixels = {
             Estimate <- as.im(NA_real_, Window(X))
             if(!univariate) {
               Estimate <- as.solist(rep(list(Estimate), nc))
               names(Estimate) <- cn
             }
             result <- if(!se) {
                         Estimate
                       } else if(univariate) {
                         solist(estimate=Estimate, SE=Estimate)
                       } else {
                         list(estimate=Estimate, SE=Estimate)
                       }
             })
    return(result)
  }

  ## ensure weights are numeric
  if(weighted <- !missing(weights) && !is.null(weights)) {
    pa <- parent.frame()
    weights <- pointweights(X, weights=weights, parent=pa)
    weighted <- !is.null(weights)
  } else weights <- NULL

  ## geometric mean smoothing
  if(geometric) 
    return(ExpSmoothLog(X, sigma=sigma, ..., at=at,
                        adjust=adjust, varcov=varcov,
                        kernel=kernel, scalekernel=scalekernel,
                        se=se, loctype=loctype, wtype=wtype,
                        weights=weights, edge=edge, diggle=diggle))

  ## determine smoothing parameters
  if(scalekernel) {
    ker <- resolve.2D.kernel(sigma=sigma, ...,
                             adjust=adjust, varcov=varcov,
                             kernel=kernel, 
                             x=X, bwfun=bw.smoothppp, allow.zero=TRUE)
    sigma <- ker$sigma
    varcov <- ker$varcov
    adjust <- 1
  } 
  
  ## ............ infinite bandwidth ...............................
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    if(!weighted) weights <- rep(1, nX)
    wtmark <- weights * marx 
    totwt <- sum(weights)
    totwtmark <- if(univariate) sum(wtmark) else colSums(wtmark)
    W <- Window(X)
    switch(at,
           pixels = {
             Estimate <- solapply(totwtmark/totwt, as.im, W=W, ...)
             names(Estimate) <- colnames(marx)
             if(univariate) Estimate <- Estimate[[1L]]
           },
           points = {
             denominator <- rep(totwt, nX)
             numerator <- rep(totwtmark, each=nX)
             if(!univariate) numerator <- matrix(numerator, nrow=nX)
             if(leaveoneout) {
               numerator <- numerator - wtmark
               denominator <- denominator - weights
             }
             Estimate <- numerator/denominator
             if(!univariate)
               colnames(Estimate) <- colnames(marx)
           })
    if(!se) {
      result <- Estimate
    } else {
      ## calculate standard error (constant value)
      if(univariate) {
        V <- if(!weighted) var(marx) else
             switch(loctype,
                    fixed = {
                      switch(wtype,
                             multiplicity = weighted.var(marx, weights),
                             importance = var(marx * weights))
                    },
                    random = {
                      switch(wtype,
                             multiplicity = VarOfWtdMean(marx, weights),
                             importance = VarOfWtdMean(marx, weights^2))
                    })
        SE <- sqrt(V) # single value
      } else {
        V <- if(!weighted) sapply(marx, var) else
             switch(loctype,
                    fixed = {
                      switch(wtype,
                             multiplicity = sapply(marx, weighted.var,
                                                   wt=weights),
                             importance = sapply(marx * weights, var))
                    },
                    random = {
                      switch(wtype,
                             multiplicity = sapply(marx, VarOfWtdMean,
                                                   weights=weights),
                             importance = sapply(marx, VarOfWtdMean,
                                                 weights=weights^2))
                    })
        SE <- sqrt(V) # vector
      }
      ## replicate constant value
      switch(at,
             pixels = {
               if(univariate) {
                 SE <- as.im(SE, W=W, ...) # constant image
                 result <- solist(estimate=Estimate, SE=SE)
               } else {
                 SE <- solapply(SE, as.im, ...) # list of images
                 result <- list(estimate=Estimate, SE=SE)
               }
             },
             points = {
               if(univariate) {
                 SE <- rep(SE, nX)
                 result <- cbind(estimate=Estimate, SE=SE)
               } else {
                 SE <- matrix(SE[col(marx)], nrow=nX, ncol=ncol(marx))
                 colnames(SE) <- colnames(marx)
                 result <- cbind(estimate=Estimate, SE=SE)
               }
             })
    }
    return(result)
  }

  ## .................. finite bandwidth .............................
  
  ## Diggle's edge correction?
  if(diggle && !edge) warning("Option diggle=TRUE overridden by edge=FALSE")
  diggle <- diggle && edge
  ##
  ## cutoff distance (beyond which the kernel value is treated as zero)
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, adjust=adjust, ...,
                           fatal=TRUE)
  
  ## ...................  bandwidth close to zero .....................
  if(cutoff < minnndist(X)) {
    # very small bandwidth
    if(!leaveoneout && at=="points") {
      warning(paste("Bandwidth is close to zero:",
                    "original values returned"))
      Estimate <- marks(X)
    } else {
      warning(paste("Bandwidth is close to zero:",
                    "nearest-neighbour interpolation performed"))
      Estimate <- nnmark(X, ..., k=1, at=at)
    }
    if(!se) {
      result <- Estimate
    } else {
      SE <- 0 * Estimate
      switch(at,
             pixels = {
               if(univariate) {
                 result <- solist(estimate=Estimate, SE=SE)
               } else {
                 result <- list(estimate=Estimate, SE=SE)
               }
             },
             points = {
               result <- cbind(estimate=Estimate, SE=SE)
             })
    }
    return(result)
  }

  ## ................... bandwidth >> 0 .........................

  if(se) {
    ## ................... STANDARD ERROR CALCULATION ..............
    ## This has to be done now because the subsequent code
    ## fiddles with the weights.
    weightspower <- if(!weighted) 1 else switch(wtype,
                                                importance = weights^2,
                                                multiplicity = weights)
    if(diggle) {
      ## Jones-Diggle correction weights e(x_i)
      edgeim <- second.moment.calc(X, sigma, what="edge", ...,
                                   varcov=varcov)
      edgeX <- safelookup(edgeim, X, warn=FALSE)
      invmassX <- 1/edgeX
      invmassX[!is.finite(invmassX)] <- 0
    } else {
      invmassX <- 1
    }
    ## 
    switch(loctype,
           random = {
             denom <- density(X, sigma=sigma, ...,
                              weights=weights,
                              edge=edge, diggle=diggle,
                              at=at, leaveoneout=leaveoneout)
             numer <- density(X, sigma=sigma, ...,
                              weights=if(weighted) weights * marx else marx,
                              edge=edge, diggle=diggle,
                              at=at, leaveoneout=leaveoneout)
             varNum <- density(X, sigma=sigma, ..., kerpow=2,
                               weights=weightspower * marx^2 * invmassX^2,
                               edge=FALSE, diggle=FALSE,
                               at=at, leaveoneout=leaveoneout)
             covND <- density(X, sigma=sigma, ..., kerpow=2,
                              weights=weightspower * marx * invmassX^2,
                              edge=FALSE, diggle=FALSE,
                              at=at, leaveoneout=leaveoneout)
             varDen <- density(X, sigma=sigma, ..., kerpow=2,
                               weights=weightspower * invmassX^2,
                               edge=FALSE, diggle=FALSE,
                               at=at, leaveoneout=leaveoneout)
             if(univariate || at == "points") {
               Vest <- DeltaMethodVarOfRatio(numer, denom,
                                          varNum, varDen, covND)
             } else {
               Vest <- mapply(DeltaMethodVarOfRatio,
                           num=numer, varnum=varNum, covnumden=covND,
                           MoreArgs = list(den=denom, varden=varDen),
                           SIMPLIFY=FALSE)
               Vest <- as.solist(Vest)
             }
           },
           fixed = {
             ## Use leave-one-out deviation
             dev <- marks(X) - Smooth(X, sigma=sigma, ...,
                                      weights=weights, edge=edge, diggle=diggle,
                                      at="points", leaveoneout=TRUE)
             if(!univariate)
               dev <- asNumericMatrix(as.data.frame(dev))
             ## calculate variance of numerator using leave-one-out estimates
             dataweight <- dev^2
             if(weighted)
               dataweight <- dataweight * switch(wtype,
                                                 importance=weights^2,
                                                 multiplicity=weights)
             if(edge && diggle) 
               dataweight <- dataweight * invmassX^2

             ## variance of numerator
             Vnum <- density(unmark(X), sigma=sigma, kerpow=2,
                             weights=dataweight,
                             at=at, leaveoneout=leaveoneout,
                             edge=FALSE, diggle=FALSE, # sic
                             positive=TRUE)
             
             if(at == "points" && !univariate)
               Vnum <- asNumericMatrix(as.data.frame(Vnum))
    
             ## rescale by denominator^2
             Den <- density(unmark(X), sigma=sigma,
                            weights=weights,
                            edge=edge && diggle, diggle=diggle, # sic
                            at=at, leaveoneout=leaveoneout, positive=TRUE)
    
             Vest <- if(at == "points") Vnum/Den^2 else
                            imagelistOp(Vnum, Den^2, "/")
                     
           })
    SE <- sqrt(Vest)
  }
  
  ##  ------------------------------------------------------------
  ##  >>>>>>>>>>>. MAIN CALCULATION OF ESTIMATE <<<<<<<<<<<<<<<<<<
  ##  ------------------------------------------------------------
  
  if(diggle) {
    ## absorb Diggle edge correction into weights vector
    edg <- second.moment.calc(X, sigma, what="edge", ...,
                              varcov=varcov, adjust=adjust,
                              kernel=kernel, scalekernel=scalekernel)
    ei <- safelookup(edg, X, warn=FALSE)
    weights <- if(weighted) weights/ei else 1/ei
    weights[!is.finite(weights)] <- 0
    weighted <- TRUE
  }
  ## rescale weights to avoid numerical gremlins
  if(weighted && ((mw <- median(abs(weights))) > 0))
    weights <- weights/mw

  ## calculate...
  uhoh <- NULL
  if(!is.data.frame(marx)) {
    # ........ vector of marks ...................
    values <- marx
    if(is.factor(values)) 
      warning("Factor valued marks were converted to integers", call.=FALSE)
    values <- as.numeric(values)
    ## detect constant values
    ra <- range(values, na.rm=TRUE)
    if(diff(ra) == 0) {
      switch(at,
             points = {
               result <- values
             },
             pixels = {
               M <- do.call.matched(as.mask, list(w=as.owin(X), ...))
               result <- as.im(ra[1], M)
             })
    } else {
      switch(at,
             points={
               result <-
                 do.call(smoothpointsEngine,
                         resolve.defaults(list(x=quote(X),
                                               values=quote(values), 
					       weights=quote(weights),
                                               leaveoneout=leaveoneout,
                                               sigma=sigma, varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
             },
             pixels={
               values.weights <- if(weighted) values * weights else values
               dont.complain.about(values.weights)
               numerator <-
                 do.call(density.ppp,
                         resolve.defaults(list(x=quote(X),
                                               at="pixels",
                                               weights = quote(values.weights),
                                               sigma=sigma, varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
               denominator <-
                 do.call(density.ppp,
                         resolve.defaults(list(x=quote(X),
                                               at="pixels",
                                               weights = quote(weights),
                                               sigma=sigma,
                                               varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
               result <- eval.im(numerator/denominator)
               ## trap small values of denominator
               ## trap NaN and +/- Inf values of result, but not NA
               eps <- .Machine$double.eps
               nbg <- eval.im(is.infinite(result)
                              | is.nan(result)
                              | (denominator < eps))
               uhoh <- attr(numerator, "warnings")
               if(any(as.matrix(nbg), na.rm=TRUE)) {
                 warning(paste("Numerical underflow detected:",
                               "sigma is probably too small"),
                         call.=FALSE)
                 uhoh <- unique(c(uhoh, "underflow"))
                 ## l'Hopital's rule
                 distX <- distmap(X, xy=numerator)
                 whichnn <- attr(distX, "index")
                 nnvalues <- eval.im(values[whichnn])
                 result[nbg] <- nnvalues[nbg]
               }
             })
    }
  } else {
    ## ......... data frame of marks ..................
    ## convert to numerical values
    if(any(sapply(as.list(marx), is.factor)))
      warning("Factor columns of marks were converted to integers", call.=FALSE)
    marx <- asNumericMatrix(marx)
    ## detect constant columns
    ra <- apply(marx, 2, range, na.rm=TRUE)
    isconst <- (apply(ra, 2, diff) == 0)
    if(anyisconst <- any(isconst)) {
      oldmarx <- marx
#      oldX <- X
      marx <- marx[, !isconst]
      X <- X %mark% marx
    }
    if(any(!isconst)) {
      ## compute denominator
      denominator <-
        do.call(density.ppp,
                resolve.defaults(list(x=quote(X),
                                      at=at,
                                      leaveoneout=leaveoneout,
                                      weights = quote(weights),
                                      sigma=sigma, varcov=varcov,
                                      kernel=kernel,
                                      scalekernel=scalekernel,
                                      edge=FALSE),
                                 list(...)))
      ## compute numerator for each column of marks
      marx.weights <- if(weighted) marx * weights else marx
      dont.complain.about(marx.weights)
      numerators <-
        do.call(density.ppp,
                resolve.defaults(list(x=quote(X),
                                      at=at,
                                      leaveoneout=leaveoneout,
                                      weights = quote(marx.weights),
                                      sigma=sigma, varcov=varcov,
                                      kernel=kernel,
                                      scalekernel=scalekernel,
                                      edge=FALSE),
                                 list(...)))
      uhoh <- attr(numerators, "warnings")
      ## calculate ratios
      switch(at,
             points={
               if(is.null(uhoh)) {
                 ## numerators is a matrix (or may have dropped to vector)
                 if(is.data.frame(numerators)) {
                   numerators <- as.matrix(numerators)
                 } else if(!is.matrix(numerators)) {
                   numerators <- matrix(unlist(numerators), nrow=npoints(X))
                 }
                 ratio <- numerators/denominator
                 if(any(badpoints <- matrowany(!is.finite(ratio)))) {
                   whichnnX <- nnwhich(X)
                   ratio[badpoints,] <-
                     as.matrix(marx[whichnnX[badpoints], , drop=FALSE])
                 }
               } else {
                 warning("returning original values")
                 ratio <- marx
               }
               result <- as.data.frame(ratio)
               colnames(result) <- colnames(marx)
             },
             pixels={
               ## numerators is a list of images (or may have dropped to 'im')
               if(is.im(numerators))
                 numerators <- list(numerators)
               result <- solapply(numerators, "/", e2=denominator)
               eps <- .Machine$double.eps
               denOK <- eval.im(denominator >= eps)
               if(!is.null(uhoh) || !all(denOK)) {
                 ## compute nearest neighbour map on same raster
                 distX <- distmap(X, xy=denominator)
                 whichnnX <- attr(distX, "index")
                 ## fix images
                 allgood <- TRUE
                 for(j in 1:length(result)) {
                   ratj <- result[[j]]
                   valj <- marx[,j]
                   goodj <- eval.im(is.finite(ratj) & denOK)
                   result[[j]] <- eval.im(goodj, ratj, valj[whichnnX])
                   allgood <- allgood && all(goodj)
                 }
                 if(!allgood) {
                   warning(paste("Numerical underflow detected:",
                                 "sigma is probably too small"),
                           call.=FALSE)
                   uhoh <- unique(c(uhoh, "underflow"))
                 }
               }
               names(result) <- colnames(marx)
             })
    } else result <- NULL 
    if(anyisconst) {
      partresult <- result
      switch(at,
             points = {
               nX <- npoints(X)
               result <- matrix(, nX, ncol(oldmarx))
               if(length(partresult) > 0)
                 result[,!isconst] <- as.matrix(partresult)
               result[,isconst] <- rep(ra[1,isconst], each=nX)
               colnames(result) <- colnames(oldmarx)
             },
             pixels = {
               result <- vector(mode="list", length=ncol(oldmarx))
               if(length(partresult) > 0) {
                 result[!isconst] <- partresult
                 M <- as.owin(partresult[[1]])
               } else {
                 M <- do.call.matched(as.mask, list(w=as.owin(X), ...))
               }
               result[isconst] <- lapply(ra[1, isconst], as.im, W=M)
               result <- as.solist(result)
               names(result) <- colnames(oldmarx)
             })
    }
  }
  ## wrap up
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  if(length(uhoh)) attr(result, "warnings") <- uhoh 
  ## tack on standard errors?
  if(se) {
    result <- list(estimate=result, SE=SE)
    switch(at,
           points = {
             result <- do.call(cbind, result)
           },
           pixels = {
             if(univariate) result <- as.solist(result)
           })
  }
  return(result)
}


smoothpointsEngine <- function(x, values, sigma, ...,
                               kernel="gaussian", 
                               scalekernel=is.character(kernel),
                               weights=NULL, varcov=NULL,
                               leaveoneout=TRUE,
                               sorted=FALSE, cutoff=NULL,
                               debug=FALSE) {
  stopifnot(is.logical(leaveoneout))

  if(!is.null(dim(values)))
    stop("Internal error: smoothpointsEngine does not support multidimensional values")
  
  #' detect constant values
  if(diff(range(values, na.rm=TRUE)) == 0) { 
    result <- values
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }

  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian")

  ## Handle weights that are meant to be null
  if(length(weights) == 0)
     weights <- NULL

  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nX <- npoints(x)
    if(is.null(weights)) weights <- rep(1, nX)
    wtval <- weights * values
    totwt <- sum(weights)
    totwtval <- sum(wtval) 
    denominator <- rep(totwt, nX)
    numerator <- rep(totwtval, nX)
    if(leaveoneout) {
      numerator <- numerator - wtval
      denominator <- denominator - weights
    }
    result <- numerator/denominator
    return(result)
  }
  
  ## cutoff distance (beyond which the kernel value is treated as zero)
  ## NB: input argument 'cutoff' is either NULL or
  ##     an absolute distance (if scalekernel=FALSE)
  ##     a number of standard deviations (if scalekernel=TRUE)
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, cutoff=cutoff,
                           fatal=TRUE)
  ## cutoff is now an absolute distance
  if(debug)
    cat(paste("cutoff=", cutoff, "\n"))
  
  # detect very small bandwidth
  nnd <- nndist(x)
  nnrange <- range(nnd)
  if(cutoff < nnrange[1]) {
    if(leaveoneout && (npoints(x) > 1)) {
      warning("Very small bandwidth; values of nearest neighbours returned")
      result <- values[nnwhich(x)]
    } else {
      warning("Very small bandwidth; original values returned")
      result <- values
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }

  if(leaveoneout) {
    # ensure cutoff includes at least one point
    cutoff <- max(1.1 * nnrange[2], cutoff)
  }

  sd <- if(is.null(varcov)) sigma else sqrt(max(eigen(varcov)$values))
  
  if(isgauss &&
     spatstat.options("densityTransform") &&
     spatstat.options("densityC")) {
    ## .................. experimental C code .....................
    if(debug)
      cat('Transforming to standard coordinates (densityTransform=TRUE).\n')
    npts <- npoints(x)
    result <- numeric(npts)
    ## transform to standard coordinates
    xx <- x$x
    yy <- x$y
    if(is.null(varcov)) {
      xx <- xx/(sqrt(2) * sigma)
      yy <- yy/(sqrt(2) * sigma)
    } else {
      Sinv <- solve(varcov)
      xy <- cbind(xx, yy) %*% matrixsqrt(Sinv/2)
      xx <- xy[,1]
      yy <- xy[,2]
      sorted <- FALSE
    }
    ## cutoff in standard coordinates
    cutoff <- cutoff/(sqrt(2) * sd)
    ## sort into increasing order of x coordinate (required by C code)
    if(!sorted) {
      oo <- fave.order(xx)
      xx <- xx[oo]
      yy <- yy[oo]
      vv <- values[oo]
    } else {
      vv <- values
    }
    if(is.null(weights)) {
      zz <- .C(SE_Gsmoopt,
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               v       = as.double(vv),
               self    = as.integer(!leaveoneout),
               rmaxi   = as.double(cutoff),
               result  = as.double(double(npts)),
               PACKAGE="spatstat.explore")
      if(sorted) result <- zz$result else result[oo] <- zz$result
    } else {
      wtsort <- weights[oo]
      zz <- .C(SE_Gwtsmoopt,
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               v       = as.double(vv),
               self    = as.integer(!leaveoneout),
               rmaxi   = as.double(cutoff),
               weight  = as.double(wtsort),
               result  = as.double(double(npts)),
               PACKAGE="spatstat.explore")
      if(sorted) result <- zz$result else result[oo] <- zz$result
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnwhich(x)[nbg]]
    }
  } else if(isgauss && spatstat.options("densityC")) {
    # .................. C code ...........................
    if(debug)
      cat('Using standard code (densityC=TRUE).\n')
    npts <- npoints(x)
    result <- numeric(npts)
    # sort into increasing order of x coordinate (required by C code)
    if(sorted) {
      xx <- x$x
      yy <- x$y
      vv <- values
    } else {
      oo <- fave.order(x$x)
      xx <- x$x[oo]
      yy <- x$y[oo]
      vv <- values[oo]
    }
    if(is.null(varcov)) {
      # isotropic kernel
      if(is.null(weights)) {
        zz <- .C(SE_smoopt,
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result  = as.double(double(npts)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C(SE_wtsmoopt,
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    } else {
      # anisotropic kernel
      Sinv <- solve(varcov)
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C(SE_asmoopt,
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(npts)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C(SE_awtsmoopt,
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnwhich(x)[nbg]]
    }
  } else {
    #' Either a non-Gaussian kernel or using older, partly interpreted code
    #' compute weighted densities
    if(debug)
      cat('Using partly-interpreted code.\n')
    if(is.null(weights)) {
      # weights are implicitly equal to 1
      numerator <- do.call(density.ppp,
                         resolve.defaults(list(x=quote(x), at="points",
                                               weights = quote(values),
                                               sigma=sigma,
                                               varcov=varcov,
                                               leaveoneout=leaveoneout,
                                               sorted=sorted,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               cutoff=cutoff),
                                          list(...),
                                          list(edge=FALSE)))
      denominator <- do.call(density.ppp,
                             resolve.defaults(list(x=quote(x), at="points",
                                                   sigma=sigma,
                                                   varcov=varcov,
                                                   leaveoneout=leaveoneout,
                                                   sorted=sorted,
                                                   kernel=kernel,
                                                   scalekernel=scalekernel,
                                                   cutoff=cutoff),
                                              list(...),
                                              list(edge=FALSE)))
    } else {
      values.weights <- values * weights
      dont.complain.about(values.weights)
      numerator <- do.call(density.ppp,
                           resolve.defaults(list(x=quote(x), at="points",
                                                 weights = quote(values.weights),
                                                 sigma=sigma,
                                                 varcov=varcov,
                                                 leaveoneout=leaveoneout,
                                                 sorted=sorted,
                                                 kernel=kernel,
                                                 scalekernel=scalekernel,
                                                 cutoff=cutoff),
                                            list(...),
                                            list(edge=FALSE)))
      denominator <- do.call(density.ppp,
                             resolve.defaults(list(x=quote(x), at="points",
                                                   weights = quote(weights),
                                                   sigma=sigma,
                                                   varcov=varcov,
                                                   leaveoneout=leaveoneout,
                                                   sorted=sorted,
                                                   kernel=kernel,
                                                   scalekernel=scalekernel,
                                                   cutoff=cutoff),
                                              list(...),
                                              list(edge=FALSE)))
    }
    if(is.null(uhoh <- attr(numerator, "warnings"))) {
      result <- numerator/denominator
      result <- ifelseXB(is.finite(result), result, NA_real_)
    } else {
      warning("returning original values")
      result <- values
      attr(result, "warnings") <- uhoh
    }
  }
  # pack up and return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


markmean <- function(X, ...) {
  stopifnot(is.marked(X))
  Y <- Smooth(X, ...)
  return(Y)
}

markvar  <- function(X, sigma=NULL, ..., weights=NULL, varcov=NULL) {
  stopifnot(is.marked(X))
  if(is.expression(weights)) 
    weights <- eval(weights, envir=as.data.frame(X), enclos=parent.frame())
  E1 <- Smooth(X, sigma=sigma, varcov=varcov, weights=weights, ...)
  X2 <- X %mark% marks(X)^2
  ## ensure smoothing bandwidth is the same!
  sigma <- attr(E1, "sigma")
  varcov <- attr(E1, "varcov")
  E2 <- Smooth(X2, sigma=sigma, varcov=varcov, weights=weights, ...)
  V <- eval.im(E2 - E1^2)
  return(V)
}

bw.smoothppp <- function(X, nh=spatstat.options("n.bandwidth"),
                         hmin=NULL, hmax=NULL, warn=TRUE,
                         kernel="gaussian", varcov1=NULL) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  if(!is.null(varcov1))
    check.nmatrix(varcov1, 2, things="spatial dimensions", mname="varcov1")
  if(is.function(kernel))
    stop("Custom kernel functions are not yet supported in bw.smoothppp")
  X <- coerce.marks.numeric(X)
  # rearrange in ascending order of x-coordinate (for C code)
  X <- X[fave.order(X$x)]
  #
  marx <- marks(X)
  dimmarx <- dim(marx)
  if(!is.null(dimmarx))
    marx <- as.matrix(as.data.frame(marx))
  # determine a range of bandwidth values
#  n <- npoints(X)
  if(is.null(hmin) || is.null(hmax)) {
    W <- Window(X)
    d <- diameter(as.rectangle(W))
    # Stoyan's rule of thumb 
    stoyan <- bw.stoyan(X)
    # rule of thumb based on nearest-neighbour distances
    nnd <- nndist(unique(X))
    if(any(ok <- is.finite(nnd) & (nnd > 0))) {
      nnd <- nnd[ok]
    } else {
      nnd <- d/16
    }
    if(!is.null(varcov1)) {
      dref <- (det(varcov1))^(1/4)
      d      <- d/dref
      stoyan <- stoyan/dref
      nnd    <- nnd/dref
    }
    if(is.null(hmin)) {
      hmin <- max(1.1 * min(nnd), stoyan/5)
      hmin <- min(d/8, hmin)
    }
    if(is.null(hmax)) {
      hmax <- max(stoyan * 20, 3 * mean(nnd), hmin * 2)
      hmax <- min(d/2, hmax)
    }
  } else stopifnot(hmin < hmax)
  #
  h <- geomseq(from=hmin, to=hmax, length.out=nh)
  cv <- numeric(nh)
  # 
  # compute cross-validation criterion
  for(i in seq_len(nh)) {
    if(is.null(varcov1)) {
      yhat <- Smooth(X, sigma = h[i],
                     at="points", leaveoneout=TRUE,
                     kernel=kernel, sorted=TRUE)
    } else {
      yhat <- Smooth(X, varcov = (h[i]^2)  * varcov1,
                     at="points", leaveoneout=TRUE,
                     kernel=kernel, sorted=TRUE)
    }
    if(!is.null(dimmarx))
      yhat <- as.matrix(as.data.frame(yhat))
    cv[i] <- mean((marx - yhat)^2)
  }

  # optimize
  result <- bw.optim(cv, h, 
                     hname="sigma",
                     creator="bw.smoothppp",
                     criterion="Least Squares Cross-Validation",
                     warnextreme=warn,
                     hargnames=c("hmin", "hmax"),
                     unitname=if(is.null(varcov1)) unitname(X) else NULL,
                     template=varcov1, exponent=2)
  return(result)
}

smoothcrossEngine <- function(Xdata, Xquery, values, sigma, ...,
                              weights=NULL, varcov=NULL,
                              kernel="gaussian", 
                              scalekernel=is.character(kernel),
                              sorted=FALSE,
                              cutoff=NULL) {

  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian") && scalekernel
  if(isTRUE(list(...)$se))
    warning("Standard errors are not yet supported", call.=FALSE)
  
  if(!is.null(dim(weights)))
    stop("weights must be a vector")

  ndata <- npoints(Xdata)
  nquery <- npoints(Xquery)
  
  if(nquery == 0 || ndata == 0) {
    if(is.null(dim(values))) return(rep(NA_real_, nquery))
    nuttin <- matrix(NA_real_, nrow=nquery, ncol=ncol(values))
    colnames(nuttin) <- colnames(values)
    return(nuttin)
  }

  # validate weights
  if(is.matrix(values) || is.data.frame(values)) {
    k <- ncol(values)
    stopifnot(nrow(values) == npoints(Xdata))
    values <- as.data.frame(values)
  } else {
    k <- 1L
    stopifnot(length(values) == npoints(Xdata) || length(values) == 1)
    if(length(values) == 1L)
      values <- rep(values, ndata)
  }

  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    if(is.null(weights)) weights <- rep(1, ndata)
    univariate <- is.null(dim(values))
    wtval <- weights * values 
    totwt <- sum(weights)
    totwtval <- if(univariate) sum(wtval) else colSums(wtval)
    denominator <- rep(totwt, nquery)
    numerator <- rep(totwtval, each=nquery)
    if(!univariate) numerator <- matrix(numerator, nrow=nquery)
    result <- numerator/denominator
    if(!univariate)
      colnames(result) <- colnames(values)
    return(result)
  }
  
  ## cutoff distance (beyond which the kernel value is treated as zero)
  ## NB: input argument 'cutoff' is either NULL or
  ##     an absolute distance (if scalekernel=FALSE)
  ##     a number of standard deviations (if scalekernel=TRUE)
  cutoff.orig <- cutoff
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, cutoff=cutoff,
                           fatal=TRUE)
  ## cutoff is now an absolute distance

  ## detect very small bandwidth
  nnc <- nncross(Xquery, Xdata)
  if(cutoff < min(nnc$dist)) {
    if(ndata > 1) {
      warning("Very small bandwidth; values of nearest neighbours returned")
      nw <- nnc$which
      result <- if(k == 1) values[nw] else values[nw,,drop=FALSE]
    } else {
      warning("Very small bandwidth; original values returned")
      result <- values
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }
  
  ## Handle weights that are meant to be null
  if(length(weights) == 0)
     weights <- NULL

  if(!isgauss) {
    ## .................. non-Gaussian kernel ........................
    close <- crosspairs(Xdata, Xquery, cutoff)
    kerij <- evaluate2Dkernel(kernel, close$dx, close$dy,
                            sigma=sigma, varcov=varcov,
                            scalekernel=scalekernel, ...)
    ## sum the (weighted) contributions
    i <- close$i # data point
    j <- close$j # query point
    jfac <- factor(j, levels=seq_len(nquery))
    wkerij <- if(is.null(weights)) kerij else kerij * weights[i]
    denominator <- tapplysum(wkerij, list(jfac))
    if(k == 1L) {
      contribij <- wkerij * values[i]
      numerator <- tapplysum(contribij, list(jfac))
      result <- numerator/denominator
    } else {
      result <- matrix(, nrow=nquery, ncol=k)
      for(kk in 1:k) {
        contribij <- wkerij * values[i, kk]
        numeratorkk <- tapplysum(contribij, list(jfac))
        result[,kk] <- numeratorkk/denominator
      }
    }
    ## trap bad values
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      ## NaN or +/-Inf can occur if bandwidth is small
      ## Use value at nearest neighbour (by l'Hopital's rule)
      nnw <- nnc$which
      if(k == 1L) {
        result[nbg] <- values[nnw[nbg]]
      } else {
        bad <- which(nbg, arr.ind=TRUE)
        badrow <- bad[,"row"]
        badcol <- bad[,"col"]
        result[nbg] <- values[cbind(nnw[badrow], badcol)]
      }
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  } 

  ## .................. Gaussian kernel henceforth ........................
  
  ## handle multiple columns of values
  if(is.matrix(values) || is.data.frame(values)) {
    k <- ncol(values)
    stopifnot(nrow(values) == npoints(Xdata))
    values <- as.data.frame(values)
    result <- matrix(, nquery, k)
    colnames(result) <- colnames(values)
    if(!sorted) {
      ood <- fave.order(Xdata$x)
      Xdata <- Xdata[ood]
      values <- values[ood, ]
      ooq <- fave.order(Xquery$x)
      Xquery <- Xquery[ooq]
    }
    for(j in 1:k) 
      result[,j] <- smoothcrossEngine(Xdata, Xquery, values[,j],
                                      sigma=sigma, varcov=varcov,
                                      weights=weights,
                                      kernel=kernel, scalekernel=scalekernel,
                                      cutoff=cutoff.orig,
                                      sorted=TRUE, ...)
    if(!sorted) {
      sortresult <- result
      result[ooq,] <- sortresult
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }

  ## values must be a vector
  stopifnot(length(values) == npoints(Xdata) || length(values) == 1)
  if(length(values) == 1) values <- rep(values, ndata)

  result <- numeric(nquery) 
  ## coordinates and values
  xq <- Xquery$x
  yq <- Xquery$y
  xd <- Xdata$x
  yd <- Xdata$y
  vd <- values
  if(!sorted) {
    ## sort into increasing order of x coordinate (required by C code)
    ooq <- fave.order(Xquery$x)
    xq <- xq[ooq]
    yq <- yq[ooq]
    ood <- fave.order(Xdata$x)
    xd <- xd[ood]
    yd <- yd[ood]
    vd <- vd[ood] 
  }
  sd <- if(is.null(varcov)) sigma else sqrt(min(eigen(varcov)$values))
  if(is.null(varcov)) {
    ## isotropic kernel
    if(is.null(weights)) {
      zz <- .C(SE_crsmoopt,
               nquery      = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               vd      = as.double(vd),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)),
               PACKAGE="spatstat.explore")
      if(sorted) result <- zz$result else result[ooq] <- zz$result
    } else {
      wtsort <- if(sorted) weights else weights[ood]
      zz <- .C(SE_wtcrsmoopt,
               nquery      = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               vd      = as.double(vd),
               wd      = as.double(wtsort),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)),
               PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      }
    } else {
      # anisotropic kernel
      Sinv <- solve(varcov)
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C(SE_acrsmoopt,
                 nquery      = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 vd      = as.double(vd),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(nquery)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      } else {
        wtsort <- if(sorted) weights else weights[ood]
        zz <- .C(SE_awtcrsmoopt,
                 nquery      = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 vd      = as.double(vd),
                 wd      = as.double(wtsort),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(nquery)),
                 PACKAGE="spatstat.explore")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      }
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnc$which[nbg]]
    }
  # pack up and return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

ExpSmoothLog <- function(X, ..., at=c("pixels", "points"), weights=NULL,
                         se=FALSE) {
  if(se)
    stop("Standard errors are not yet supported when geometric=TRUE")
  verifyclass(X, "ppp")
  at <- match.arg(at)
  if(!is.null(weights)) 
    check.nvector(weights, npoints(X), vname="weights")
  X <- coerce.marks.numeric(X)
  marx <- marks(X)
  d <- dim(marx)
  if(!is.null(d) && d[2] > 1) {
    switch(at,
           points = {
             Z <- lapply(unstack(X), ExpSmoothLog, ...,
                         at=at, weights=weights)
             Z <- do.call(data.frame, Z)
           },
           pixels = {
             Z <- solapply(unstack(X), ExpSmoothLog, ...,
                           at=at, weights=weights)
           })
    return(Z)
  }
  # vector or single column of numeric marks
  v <- as.numeric(marx)
  vmin <- min(v)
  if(vmin < 0) stop("Negative values in geometric mean smoothing",
                       call.=FALSE)
  Y <- X %mark% log(v)
  if(vmin > 0) {
    Z <- Smooth(Y, ..., at=at, weights=weights)
  } else {
    yok <- is.finite(marks(Y))
    YOK <- Y[yok]
    weightsOK <- if(is.null(weights)) NULL else weights[yok]
    switch(at,
           points = {
             Z <- rep(-Inf, npoints(X))
             Z[yok] <- Smooth(YOK, ..., at=at, weights=weightsOK)
           },
           pixels = {
             isfinite <- nnmark(Y %mark% yok, ...)
             support <- solutionset(isfinite)
             Window(YOK) <- support
             Z <- as.im(-Inf, W=Window(Y), ...)
             Z[support] <- Smooth(YOK, ..., at=at, weights=weightsOK)[]
           })
  }
  return(exp(Z))
}

VarOfWtdMean <- function(marx, weights) {
  ## weighted average
  totwt     <- sum(weights)
  totwtmark <- sum(weights * marx)
  Estimate  <- totwtmark/totwt
  ## delta method approximation to variance of weighted average
  varnum    <- sum(weights * marx^2)
  varden    <- totwt
  covnumden <- totwtmark
  V <- DeltaMethodVarOfRatio(totwtmark, totwt, varnum, varden, covnumden)
  return(V)
}

DeltaMethodVarOfRatio <- function(num, den, varnum, varden, covnumden) {
  Estimate <- num/den
  V <- Estimate^2 * (
    varnum/num^2
    - 2 * covnumden/(num * den)
    + varden/den^2
  )
  return(V)
}
