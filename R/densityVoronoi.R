#'
#'  densityVoronoi.R
#'
#'  $Revision: 1.32 $   $Date: 2026/08/13 06:07:17 $
#'


densityVoronoi <- function(X, ...) {
  UseMethod("densityVoronoi")
}

densityVoronoi.ppp <- function(X, f=1, ...,
                               counting=FALSE,
                               fixed=FALSE,
                               nrep=1, verbose=TRUE, internal=NULL) {
  stopifnot(is.ppp(X))
  if(is.NAobject(X)) return(NAobject("im"))
  X <- unmark(X)
  nX <- npoints(X)
  check.1.real(f)
  if(badprobability(f))
    stop("f should be a probability between 0 and 1")
  check.1.integer(nrep)
  stopifnot(nrep >= 1)
  anydupes <- internal$anydupes %orifnull% anyDuplicated(X, rule="deldir")
  dupes <-    internal$dupes %orifnull% (if(anydupes) duplicated(X, rule="deldir") else NULL)
  W <- Window(X)
  areaW <- area(W)
  #' ------- handle trivial cases where replication has no effect ----------
  if(nX == 0) {
    #' empty dataset
    return(as.im(0, W, ...))
  } else if(f == 0 || (fixed && floor(f * nX) == 0)) {
    ## subset of points used for tessellation is always empty
    ## Fall back on constant intensity estimate
    lam <- nX/areaW
    return(as.im(lam, W, ...))
  } else if(f == 1) {
    ## Voronoi/Dirichlet estimate
    if(!anydupes) {
      tes <- dirichlet(X)
      tesim <- nnmap(X, what="which", ...)
      num <- 1
    } else {
      UX <- X[!dupes]
      tes <- dirichlet(UX)
      tesim <- nnmap(UX, what="which", ...)
      idx <- nncross(X, UX, what="which")
      num <- as.integer(table(factor(idx, levels=seq_len(npoints(UX)))))
    }
    lam <- num/tile.areas(tes)
    out <- eval.im(lam[tesim])
    if(anyNA(out)) warning("NA pixel values generated", call.=FALSE)
    return(out)
  }
  #' ----------------- general case -----------------------------
  if(nrep > 1) {
    ## estimate is the average of nrep randomised estimates
    total <- 0
    if(verbose)
      cat(paste("Computing", nrep, "intensity estimates..."))
    state <- list()
    dupinfo <- list(anydupes=anydupes, dupes=dupes)
    for(i in seq_len(nrep)) {
      estimate <- densityVoronoi.ppp(X, f, ...,
                                     counting=counting, fixed=fixed, nrep=1,
                                     internal=dupinfo)
      total <- eval.im(total + estimate)
      if(verbose) state <- progressreport(i, nrep, state=state)
    }
    if(verbose) cat("Done.\n")
    average <- eval.im(total/nrep)
    return(average)
  }
  ## ------ This is the main calculation -------
  ## perform thinning 
  if(fixed) {
    #' random sample of fixed size 'ntess' (not equal to 0 or nX)
    ntess <- floor(f * nX)
    itess <- sample(seq_len(nX), ntess, replace=FALSE)
    nXT <- ntess
    tessfrac <- as.numeric(ntess)/nX
  } else {
    #' random thinning, retention probability f
    itess <- seq_len(nX)[thinjump(nX, f)]
    nXT <- length(itess)
    tessfrac <- f
  }
  nCount <- if(counting) nX - nXT else NULL
  #' number of points forming the tessellation
  nT <- nXT
  if(nT > 0) {
    #' form the subsample and remove duplicate points for tessellation
    Xtess <- X[itess]
    anydupesNow <- FALSE
    if(anydupes) {
      XtessU <- unique(Xtess, rule="deldir")
      nTU <- npoints(XtessU)
      anydupesNow <- (nTU < nT)
      if(anydupesNow) {
        XtessDup <- Xtess
        Xtess    <- XtessU
        nT <- nTU
      }
    }
  }
  ## determine estimation rule
  if(!counting) {
    #' Voronoi estimator of intensity of *subsample*
    totaltilemass <- nXT
    #' multiply by 1/(sampling fraction) to get intensity of X
    expansion <- 1/tessfrac
  } else {
    #' counting estimator of intensity of *un-sampled points*
    totaltilemass <- nCount
    #' multiply by 1/(*complement* of sampling fraction) to get intensity of X
    expansion <- 1/(1-tessfrac)
  }
  ## ----------- compute estimate ----------------------------
  validrule <- (nT > 0) && (totaltilemass > 0) && is.finite(expansion)
  if(!validrule) {
    ## rule collapses (tessellation undefined, divide-by-zero, etc)
    ## fall back to uniform estimate
    lam <- nX/areaW
    result <- as.im(lam, W, ...)
  } else if(nXT == 1) {
    ## rule is valid
    ## tessellation is a single tile - handle separately for efficiency
    ## spatially constant intensity estimate for X
    lam <- expansion * totaltilemass/areaW
    result <- as.im(lam, W, ...)
  } else {
    ## rule is valid
    ## general case
    tes <- dirichlet(Xtess)
    ## determine un-expanded mass in each tile
    if(!counting) {
      #' Voronoi estimator of intensity of *subsample*
      if(!anydupesNow) {
        tilemass <- 1 # each tile has 1 point
      } else {
        ## each tile may contain several points at the same location
        idx <- nncross(XtessDup, Xtess, what="which")
        tilemass <- as.integer(table(factor(idx, levels=seq_len(nT))))
      }
    } else {
      #' counting estimator of intensity of *un-sampled points*
      Xcount <- X[-itess]
      tilemap <- tileindex(Xcount$x, Xcount$y, tes)
      tilemass <- as.numeric(table(tilemap))
    }
    ## expand to give final intensity estimate for X in each tile
    lami <- expansion * tilemass/tile.areas(tes)
    ## estimate of intensity of X at each location
    tesim <- nnmap(Xtess, what="which", ...)
    result <- eval.im(lami[tesim])
  }
  if(anyNA(result)) warning("NA pixel values generated", call.=FALSE)
  return(result)
}
