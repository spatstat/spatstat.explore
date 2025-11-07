##
##
##  roc.R
##
##  Calculate ROC curve
##       (in spatstat.explore)
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner/Suman Rakshit
##
##

roc <- function(X, ...) { UseMethod("roc") }

  roc.ppp <-
  function(X, covariate, 
                               ...,
                               baseline=NULL,
                               high = TRUE,
                               weights = NULL,
                               observations = c("exact", "presence"),
                               method = "raw",
                               CI = "none", alpha=0.05,
                               subset=NULL) {
  callframe <- parent.frame()
  stopifnot(is.ppp(X))
  observations <- match.arg(observations)
  ## estimation method
  method <- match.arg(method, c("raw", "monotonic", "smooth", "rhohat", "all"),
                      several.ok = TRUE)
  method <- sub("rhohat", "smooth", method)
  if("all" %in% method)
    method <- c("raw", "monotonic", "smooth", "all")
  ## Resolve null model / baseline
  baseline.is.pattern <- is.ppp(baseline) || is.lpp(baseline)
  if(!baseline.is.pattern) {
    ## baseline is NULL, or a function/image/...  that determines null model
    nullmodel <- resolveNullModel(baseline, X, observations, ...)
  } else {
    ## baseline is a set of dummy points
    if(observations == "presence") {
      ## discretise both patterns
      X <- do.call.matched(discretise, 
                           list(X=quote(X), ..., move.points=TRUE),
                           envir=callframe)
      baseline <-
        do.call.matched(discretise,
                        list(X=quote(baseline), ..., move.points=TRUE),
                        envir=callframe)
    }
  }
  ## confidence intervals using which estimate?
  CI <- match.arg(CI, c("none", "raw", "monotonic", "smooth", "rhohat"))
  CI <- sub("rhohat", "smooth", CI)
  if(CI != "none") {
    method <- union(method, CI)
    if(baseline.is.pattern && CI != "smooth")
      warning(paste("Confidence interval for", sQuote(CI), "method",
                    "treats the baseline ('control') point pattern as fixed"),
              call.=FALSE)
  }
  #' Get covariates
  covariate <- digestCovariates(covariate, W = Window(X))
  #' ensure 'high' is a vector
  ncov <- length(covariate)
  if(length(high) == 1)
    high <- rep(high, ncov)
  stopifnot(length(high) == ncov)
  #' process
  result <- vector(mode="list", length=ncov)
  for(i in seq_len(ncov)) {
    if(baseline.is.pattern){
      result[[i]] <- rocDummy(X, baseline, covariate[[i]],
                              method=method,
                              high = high[i],
                              weights=weights,
                              subset=subset,
                              CI=CI, alpha=alpha,
                              ...)
    } else{
      result[[i]] <- rocEngine(covariate[[i]],
                               nullmodel,
                               method = method,
                               high=high[i],
                               weights=weights,
                               subset=subset,
                               CI=CI, alpha=alpha,
                               ...)
      }
  }
  names(result) <- names(covariate)
  result <- as.anylist(result)
  if(ncov == 1)
    result <- result[[1]]
  return(result)
}

roc.im <- function(X, covariate, ..., high=TRUE) {
  rocIm(X, covariate, ..., high=high)
}

rocIm <- function(X, covariate, ..., high=TRUE,
                  p=seq(0,1,length=plength),
                  plength=1024) {
  ignored <- c("method", "CI")
  if(any(found <- !is.na(match(ignored, names(list(...)))))) 
    warning(paste(ngettext(sum(found), "Argument", "Arguments"),
                  commasep(sQuote(ignored[found])),
                  ngettext(sum(found), "was", "were"),
                  "ignored by roc.im"),
            call.=FALSE)
  Z <- covariate
  FZ <- spatialcdf(Z, ..., normalise=TRUE)
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
    plength <- length(p)
  }
  if(high) p <- 1-p  # reverse again at the end
  thresh <- quantile(FZ, p)
  ## make unique threshold breaks for later use
  if(anydup <- anyDuplicated(thresh)) {
    good <- !duplicated(thresh)
    uthresh <- thresh[good]
  } else {
    uthresh <- thresh
  }
  #' ensure Z is now an image
  if(!is.im(Z)) {
    W <- Window(X)
    Z <- digestCovariates(Z, W=W)[[1L]]
    Z <- as.im(Z, W=W)
  }
  #' extract pixel values on common raster
  dat <- pairs.im(X, Z, plot=FALSE)
  xval <- dat[,"X"]
  zval <- dat[,"Z"]
  #' accumulate using 'cut' with unique breaks
  zgrp <- cut(zval,
              uthresh,
              labels=1:(length(uthresh)-1),
              include.lowest=TRUE)
  xmass <- as.numeric(tapplysum(xval, list(zgrp)))
  xcum <- if(high) c(revcumsum(xmass), 0) else c(0, cumsum(xmass))
  xcum <- xcum/max(xcum)
  if(anydup) {
    uxcum <- xcum
    xcum <- numeric(plength)
    xcum[good] <- uxcum
    if(high) {
      # ensure nonincreasing
      xcum[-plength] <- pmax(xcum[-1], xcum[-plength])
    } else {
      # ensure nondecreasing
      xcum[-1] <- pmax(xcum[-1], xcum[-plength])
    }
  }
  df <- data.frame(p = p, thresh=thresh, null=p, raw=xcum)
  if(high) df <- df[plength:1, ]
  as.roc.data.frame(df)
}

#' Temporary code provides soon-to-be-deprecated functions rocData, rocModel

rocData <- function(covariate, nullmodel, ..., high=TRUE,
                    p=seq(0, 1, length=1024)) {
  rocEngine(covariate, nullmodel, high=high, p=p, covtype="covariate")
}

rocModel <- function(lambda, nullmodel, ..., high=TRUE,
                     p=seq(0, 1, length=1024),
                     lambdatype=c("intensity", "probabilities")) {
  lambdatype <- match.arg(lambdatype)
  covtype <- if(lambdatype == "probabilities") "probability" else "intensity"
  if(!high)
    warning("Argument 'high' is ignored for point process models",
            call.=FALSE)
  rocEngine(lambda, nullmodel, high=TRUE, p=p, covtype=covtype)
}



resolveNullModel <- function(baseline, X,
                             observations=c("exact", "presence"), ...) {

  observations <- match.arg(observations) # used only when X is a point pattern
  
  if(is.ppp(baseline) || is.lpp(baseline))
    return(NULL)
  if(inherits(baseline,
              c("ppm", "kppm", "slrm", "dppm", "lppm",
                "exactppm", "exactlppm")))
    return(baseline)

  Xispattern <- is.ppp(X) || is.lpp(X)
  Xismodel <- inherits(X, c("ppm", "kppm", "slrm", "dppm", "lppm"))
  if(!(Xispattern || Xismodel))
    stop("Internal error: format of 'X' is not understood")

  if(is.null(baseline)) {
    ## baseline was not given; defaults to CSR
    if(Xispattern) 
      switch(observations,
             exact = {
               if(is.ppp(X)) return(exactppm(X))
               if(is.lpp(X)) {
                 needpackage("spatstat.linnet")
                 return(spatstat.linnet::exactlppm(X))
               }
             },
             presence = {
               if(is.ppp(X)) {
                 needpackage("spatstat.model")
                 return(spatstat.model::slrm(X ~ 1, ...))
               }
               if(is.lpp(X))
                 stop("Presence/absence analysis is not yet supported on networks")
             })
    ## X must be a model
    ensureModelSupport(X)
    return(update(X, . ~ 1))
  }

  if(is.im(baseline) || is.function(baseline) || is.owin(baseline) ||
     is.tess(baseline) || is.numeric(baseline)) {
    if(Xispattern)
      switch(observations,
             exact = {
               ## WAS: if(is.ppp(X)) return(ppm(X ~ offset(log(baseline))))
               if(is.ppp(X)) return(exactppm(X, baseline=baseline))
               if(is.lpp(X)) {
                 needpackage("spatstat.linnet")
                 return(spatstat.linnet::lppm(X ~ offset(log(baseline))))
               }
             },
             presence = {
               if(is.ppp(X)) {
                 needpackage("spatstat.model")
                 return(spatstat.model::slrm(X ~ offset(log(baseline)), ...))
               }
               if(is.lpp(X))
                 stop("Presence/absence analysis is not yet supported on networks")
             })
    ## X must be a model
    ensureModelSupport(X)
    ## expand the environment to include 'baseline'
    envy <- list2env(list(baseline=baseline),
                     parent=environment(terms(X)))
    if(inherits(X, "slrm")) {
      ## this only seems to work for 'slrm' objects 
      return(update(X, . ~ offset(log(baseline)), env=envy))
    } else {
      ## workaround for ppm and lppm
      cl <- getCall(X)
      cname <- as.character(cl[[1]])
      if(!(cname %in% c("ppm.formula", "lppm.formula")))
        stop(paste("Sorry, models fitted by", dQuote(cname),
                   "cannot be handled yet;",
                   "please refit using an explicit formula"))
      foname <- if(inherits(X, "ppm")) "Q" else "X"
      oldfo <- cl[[foname]]
      newfo <- update(as.formula(oldfo), . ~ offset(log(baseline)))
      environment(newfo) <- envy
      cl[[foname]] <- newfo
      result <- eval(cl, envir=envy)
      return(result)
    }
  }
  stop("Format of 'baseline' is not understood")
}

rocEngine <- function(discrim, nullmodel,
                      ...,
                      covtype = c("covariate", "intensity", "probability"),
                      fittedmodel = NULL,
                      method = "raw", high = TRUE, weights=NULL,
                      discrimAtPoints = NULL,
                      p = seq(0, 1, length=plength),
                      plength = 1024,
                      interpolate=FALSE, jitter=FALSE,
                      subset=NULL, bw="nrd0", adjust=1,
                      CI="none", alpha=0.05, degfreefun=Inf,
                      leftoneout=FALSE) {
  
  ##  `discrim` is the discriminant function
  ##      (an image, or a list of images for each type of point)
  ##  `fittedmodel' is NULL, or a fitted model
  ##  'discrimAtPoints' is an optional vector, providing different values for
  ##      the discriminant at the data point locations.

  #' validate arguments
  covtype <- match.arg(covtype)
  method <- match.arg(method, c("raw", "monotonic", "smooth", "rhohat", "all"),
                      several.ok = TRUE)
  method <- sub("rhohat", "smooth", method)
  if("all" %in% method)
    method <- c("raw", "monotonic", "smooth", "all")

  #' confidence interval (for one of the estimates)
  CI <- match.arg(CI, c("none", "raw", "monotonic", "smooth", "rhohat"))
  CI <- sub("rhohat", "smooth", CI)
  if(CI != "none")
    method <- union(method, CI)

  ensureModelSupport(nullmodel)

  ## Extract data pattern
  X <- spatstat.model::response(nullmodel)
  nX <- npoints(X)

  if(inherits(nullmodel, "slrm") && is.null(nullmodel$Data$dataAtPoints))
    X <- discretise(X, xy=nullmodel$Data$W, move.points=TRUE)

  ## Validate weights 
  if(!is.null(weights)) {
    if(length(weights) == 1) weights <- rep(weights, npoints(X))
    stopifnot(length(weights) == npoints(X))
    stopifnot(all(weights >= 0))
    if(!is.null(subset)) 
      weights <- weights[ppsubset(X, subset)]
    #' normalise weights to sum to 1
    weights <- weights/sum(weights)
  }

  ## Compute basic needed quantities
  d <- spatialCDFframe(nullmodel, discrim, ...,
                       covariateAtPoints = discrimAtPoints,
                       interpolate=interpolate, jitter=jitter,
                       subset=subset,
                       raster.action="ignore",
                       make.quantile.function=TRUE)
  dval <- d$values
  U <- dval$U
  UU <- if(high) 1 - U else U
  ec <- ewcdf(UU, weights=weights)
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
  }
  FZ <- dval$FZ
  FZinverse <- dval$FZinverse
  ## thresholds
  thresh <- if(high) FZinverse(1-p) else FZinverse(p)
  ## Initialize output object
  df <- data.frame(p=p, thresh=thresh, null=p)

  # Add raw estimate if requested
  if("raw" %in% method){
    df$raw <- Estimate <- ec(p)
    if(CI == "raw") {
      se <- sqrt(Estimate * (1-Estimate)/nX)
      crit <- stats::qnorm(1-alpha/2)
      df <- cbind(df,
                  data.frame(se = se,
                             lo = pmax(0.0, Estimate - crit * se),
                             hi = pmin(1.0, Estimate + crit * se)))
      method <- c(method, "se", "lo", "hi")
    }
  }

  if(!is.null(fittedmodel) || covtype != "covariate") {
    ## Add 'theoretical' curve predicted by model
    switch(covtype,
           probability = {
             ## traditional usage
             ## Discriminant is predicted presence probability p_j
             discrimvalues <- dval$Zvalues
             ## Pixel weight = presence probability
             pixelweights <- dval$EdN
           },
           intensity = {
             ## continuous version of traditional usage
             ## Discriminant is predicted intensity lambda(u)
             discrimvalues <- dval$Zvalues
             ## Pixel weight = lambda(u) du 
             pixelweights <- dval$EdN
           },
           covariate = {
             ## Non-traditional usage: model-predicted ROC of another covariate.
             ## Obtain values of fitted intensity and discriminant at each pixel
             b <- spatialCDFframe(fittedmodel, discrim, ...,
                                  interpolate=interpolate, jitter=jitter,
                                  subset=subset,
                                  raster.action="ignore")
             vb <- b$values
             discrimvalues <- vb$Zvalues
             ## Pixel weight = lambda(u) du or probability
             pixelweights  <- vb$EdN
           })
    if(high) {
      F1negZ <- ewcdf(-discrimvalues, pixelweights/sum(pixelweights))
      df$theo <- F1negZ(-FZinverse(1-p))
    } else {
      F1Z <- ewcdf(discrimvalues, pixelweights/sum(pixelweights))
      df$theo <- F1Z(FZinverse(p))
    }
  }

  ## Add smooth estimate if requested
  if("smooth" %in% method){
    doCI <- (CI == "smooth")
    est <- rocSmoothCalc(ZX       = dval$ZX,
                         ZU       = dval$Zvalues,
                         weightsX = weights,
                         weightsU = dval$EdN,
                         high=high, p=p,
                         bw=bw, adjust=adjust,
                         doCI=doCI, alpha=alpha,
                         degfreeU=degfreefun)
    if(doCI) {
      df <- cbind(df, est[,c("smooth", "se", "lo", "hi")])
      method <- c(method, "se", "lo", "hi")
    } else {
      df$smooth <- est$smooth
    }
  }

  ## Add monotonic estimate if requested
  if("monotonic" %in% method){
    f <- monotonicRhoFun(X, discrim, increasing = high,
                   subset=subset, weights=weights, baseline=nullmodel)
    est <- roc.function(f, high = high, tp = data.frame(p, thresh),
                        method = "monotonic")
    df$monotonic <- Estimate <- est$monotonic
    if(CI == "monotonic") {
      se <- sqrt(Estimate * (1-Estimate)/nX)
      crit <- stats::qnorm(1-alpha/2)
      df <- cbind(df,
                  data.frame(se = se,
                             lo = pmax(0.0, Estimate - crit * se),
                             hi = pmin(1.0, Estimate + crit * se)))
      method <- c(method, "se", "lo", "hi")
    }
  }

  ## Convert to roc+fv object with nice labels etc.
  rslt <- as.roc.data.frame(df, method = method, covtype=covtype, CI=CI,
                            leftoneout=leftoneout)

  return(rslt)
}

as.roc.data.frame <- function(x, method = NULL,
                              covtype=c("covariate", "intensity",
                                        "probability"),
                              CI="none", leftoneout=FALSE){
  ## Convert table of calculated ROC values
  ## to an object of class 'roc', 'fv' with nice labels etc.
  covtype <- match.arg(covtype)
  nam <- names(x)
  if(!all(c("p", "thresh") %in% nam))
    stop("x must have columns named 'p' and 'thresh'")
  knownestimators <- c("raw", "monotonic", "smooth", "function", "fun")
  knownextras     <- c("se", "lo", "hi")
  knownmethods  <- c(knownestimators, knownextras)
  unchangedTags <- c("p", "null", "theo", "thresh", "fun")
  if(is.null(method) || "all" %in% method) {
    ## infer method from column names
    method <- intersect(knownmethods, nam)
    ## map abbreviations
    method[method == "function"] <- "fun"
  } else {
    ## validate method(s)
    h <- match(method, knownmethods)
    if(any(unknown <- is.na(h))) {
      nuk <- sum(unknown)
      stop(paste("Unrecognised",
                 ngettext(nuk, "method", "methods"),
                 sQuote(commasep(method[unknown]))),
           call.=FALSE)
    }
    ## map abbreviations
    method[method == "function"] <- "fun"
    ## ensure columns required for method are present in data
    if(any(unavail <- !(method %in% nam))) {
      nun <- sum(unavail)
      stop(paste("Cannot perform",
                 ngettext(nun, "method", "methods"),
                 sQuote(commasep(method[unavail])),
                 "because the data were not provided in x"),
           call.=FALSE)
    }
  }
  #' sort methods in order given by 'knownmethods'
  method <- knownmethods[sort(match(method, knownmethods))]
  ## Initialize output object in correct order
  df <- data.frame(p=x$p, thresh=x$thresh, null=x$p)
  threshname <- paste("threshold for",
                switch(covtype, probability = "presence probability", covtype))
  desc <- c("fraction of area",
            threshname,
            "expected fraction of points if no effect")
  labl <- c("p", "threshold", "%s[null](p)")
  dotnames <- "null"

  ## Add theoretical curve in case of a fitted model
  if("theo" %in% nam){
    df$theo <- x$theo
    desc <- c(desc, "expected fraction of points according to model")
    labl <- c(labl, makeRocFlabel("theo"))
    dotnames <- c("theo", dotnames)
  }

  ## Add raw estimate if requested
  if("raw" %in% method){
    df[[makeRocTag("raw", leftoneout)]] <- x$raw
    desc <- c(desc,
              makeRocDesc("observed fraction of points (raw estimate)",
                       leftoneout))
    labl <- c(labl, makeRocFlabel("raw", leftoneout))
  }

  ## Add function estimate if requested
  if("fun" %in% method){
    df$fun <- x$fun
    desc <- c(desc, "observed fraction of points (function estimate)")
    labl <- c(labl, makeRocFlabel("fun"))
  }

  ## Add monotonic estimate if requested
  if("monotonic" %in% method){
    df[[makeRocTag("monotonic", leftoneout)]] <- x$monotonic
    desc <- c(desc,
              makeRocDesc("observed fraction of points (monotonic estimate)",
                       leftoneout))
    labl <- c(labl, makeRocFlabel("monotonic", leftoneout))
  }

  ## Add smooth estimate if requested
  if("smooth" %in% method){
    df[[makeRocTag("smooth", leftoneout)]] <- x$smooth
    desc <- c(desc,
              makeRocDesc("observed fraction of points (smooth estimate)",
                       leftoneout))
    labl <- c(labl, makeRocFlabel("smooth", leftoneout))
  }
  if("se" %in% method){
    df[[makeRocTag("se", leftoneout)]] <- x$se
    desc <- c(desc,
              makeRocDesc(paste("standard error of", CI, "estimate"),
                       leftoneout))
    labl <- c(labl, paste0("se", paren(makeRocFlabel(CI, leftoneout))))
  }
  if("lo" %in% method){
    df[[makeRocTag("lo", leftoneout)]] <- x$lo
    desc <- c(desc,
              makeRocDesc("lower limit of 95%% confidence band", leftoneout))
    labl <- c(labl, makeRocFlabel("lo", leftoneout))
  }
  if("hi" %in% method){
    df[[makeRocTag("hi", leftoneout)]] <- x$hi
    desc <- c(desc,
              makeRocDesc("upper limit of 95%% confidence band", leftoneout))
    labl <- c(labl, makeRocFlabel("hi", leftoneout))
  }


  ## Ensure data are sorted in increasing order of p
  nr <- nrow(df)
  if(with(df, p[1] > p[nr]))
    df <- df[nr:1, ]
  
  result <- fv(df,
               argu  = "p",
               ylab  = quote(roc(p)),
               valu  = makeRocTag(method[1], leftoneout),
               fmla  = . ~ p,
               desc  = desc,
               labl  = labl,
               fname = "roc")
  
  #' add plot info
  dotnames <- c(method, dotnames)
  if(all(c("lo", "hi") %in% method)) {
    shadenames <- c("lo", "hi")
    dotnames <- setdiff(dotnames, c("se", "lo", "hi"))
  } else {
    shadenames <- NULL
  }
  ## 
  changing <- !(dotnames %in% unchangedTags)
  dotnames[changing] <- sapply(dotnames[changing],
                               makeRocTag, leftoneout=leftoneout)
  fvnames(result, ".") <- as.character(dotnames)
  if(!is.null(shadenames)) {
    shadenames <- sapply(shadenames, makeRocTag, leftoneout=leftoneout)
    fvnames(result, ".s") <- as.character(shadenames)
  }
  #' 
  class(result) <- unique(c("roc", class(result)))
  return(result)
}

plot.roc <- function(x, fmla, ..., main, threshold=FALSE) {
  if(missing(main)) main <- short.deparse(substitute(x))
  if(missing(fmla)) fmla <- NULL
  result <- plot.fv(x, fmla=fmla, main=main, ...)
  if(is.null(fmla) && threshold) {
    p <- with(x, .x)
    thresh <- with(x, thresh)
    fobs <- with(x, .y)
    FZXinv <- approxfun(thresh, fobs, rule=2)
    FZinv <- approxfun(thresh, p, rule=2)
    tval <- prettyweird(thresh, fobs)$x
    pval <- FZXinv(tval)
    axis(4, at=pval, labels=tval)
    mtext("Threshold value", 4, 2)
    tval <- prettyweird(thresh, p)$x
    pval <- FZinv(tval)
    nmain <- sum(nzchar(main))
    nline <- if(nmain == 0) 0 else (nmain + 2)
    axis(3, at=pval, labels=tval, line = nline)
    mtext("Threshold value", 3, nline+2)
  }
  return(invisible(result))
}

rocDummy <- function(X, U, covariate, ...,
                     high=TRUE,
                     method = c("raw", "monotonic", "smooth"),
                     subset=NULL, weights=NULL, weightsU=NULL,
                     p = seq(0, 1, length=plength),
                     plength = 1024,
                     bw="nrd0", adjust=1,
                     CI="none", alpha=0.05, degfreefun=NULL) {
  method <- match.arg(method, several.ok=TRUE)
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
  }
  #' confidence interval (for one of the estimates)
  CI <- match.arg(CI, c("none", "raw", "monotonic", "smooth"))
  if(CI != "none") method <- union(method, CI)
  #' X and U are data and dummy point patterns, respectively
  if(!is.null(weights)) {
    if(length(weights) == 1) weights <- rep(weights, npoints(X))
    stopifnot(length(weights) == npoints(X))
    stopifnot(all(weights >= 0))
  }
  if(!is.null(weightsU)) {
    if(length(weightsU) == 1) weightsU <- rep(weights, npoints(U))
    stopifnot(length(weightsU) == npoints(U))
    stopifnot(all(weightsU >= 0))
  }
  if(!is.null(subset)) {
    if(!is.null(weights)) weights <- weights[ppsubset(X, subset)]
    if(!is.null(weightsU)) weightsU <- weightsU[ppsubset(U, subset)]
    X <- X[subset]
    U <- U[subset]
  }
  nX <- npoints(X)
  if(!is.null(weights)) weights <- weights/sum(weights)
  if(!is.null(weightsU)) weightsU <- weightsU/sum(weightsU)
  ZX <- evaluateCovariate(covariate, X)
  ZU <- evaluateCovariate(covariate, U)
  FZX <- ewcdf(ZX, weights=weights)
  FZU <- ewcdf(ZU, weights=weightsU)
  efzu <- environment(FZU)
  FZUinverse <- approxfun(get("y", envir=efzu),
                          get("x", envir=efzu),
                          rule=2)
  thresh <- FZUinverse(if(high) (1-p) else p)
  df <- data.frame(p=p, thresh=thresh, null=p)
  if("raw" %in% method) {
    raw <- if(high) 1 - FZX(thresh) else FZX(thresh)
    df$raw <- raw
    if(CI == "raw") {
      se <- sqrt(raw * (1-raw)/nX)
      crit <- stats::qnorm(1-alpha/2)
      df <- cbind(df,
                  data.frame(se = se,
                             lo = pmax(0.0, raw - crit * se),
                             hi = pmin(1.0, raw + crit * se)))
      method <- c(method, "se", "lo", "hi")
    }
  }
  if("smooth" %in% method) {
    doCI <- (CI == "smooth")
    est <- rocSmoothCalc(ZX, ZU, weightsX=weights, weightsU=weightsU,
                         p=p, high=high, bw=bw, adjust=adjust,
                         doCI=doCI, alpha=alpha, degfreeU=degfreefun)
    if(doCI) {
      df <- cbind(df, est[,c("smooth", "se", "lo", "hi")])
      method <- c(method, "lo", "hi")
    } else {
      df$smooth <- est$smooth
    }
  }
  if("monotonic" %in% method){
    f <- monotonicRhoFunCalc(x=ZX, massx=weights, z=ZU, weightz=weightsU,
                       increasing = high)
    est <- roc.function(f, high = high, tp = df[,c("p", "thresh")],
                        method = "monotonic")
    df$monotonic <- Estimate <- est$monotonic
    if(CI == "monotonic") {
      se <- sqrt(Estimate * (1-Estimate)/nX)
      crit <- stats::qnorm(1-alpha/2)
      df <- cbind(df,
                  data.frame(se = se,
                             lo = pmax(0.0, Estimate - crit * se),
                             hi = pmin(1.0, Estimate + crit * se)))
      method <- c(method, "se", "lo", "hi")
    }
  }
  return(as.roc.data.frame(df, method, CI=CI))
}

makeRocFlabel <- function(sub="raw", leftoneout=FALSE,
                           super=NULL, f="%s", argu="p") {
  if(is.null(super) && leftoneout) super <- "bold(\"-\")"
  z <- if(is.null(super)) {
         paste0(f, "[", sub, "](", argu, ")")
       } else {
         paste0("{", f, "[", sub, "]^{", super, "}}(", argu, ")")
       }
  return(z)
}

makeRocTag <- function(sub="raw", leftoneout=FALSE) {
  paste0(sub, if(leftoneout) "loo" else "")
}

makeRocDesc <- function(desc, leftoneout=FALSE) {
  paste0(desc, if(leftoneout) " (leave-one-out)" else "")
}




#' Compute the expected ROC curve for covariate Z
#' given that the true intensity is lambda(u) = rho(Z(u))

roc.rhohat <- function(X, ..., high = TRUE){
  covar <- attr(X, "stuff")$Zimage
  f <- as.function(X, extrapolate = TRUE)
  roc.function(f, covar, ..., high = high, method = "smooth")
}

roc.function <- function(X, covariate, ...,
                         high = TRUE, tp=NULL, method = "function",
                         nsteps = 1024){
  if(is.null(tp)) {
    if(!is.im(covariate)) covariate <- as.im(covariate, ...)
    if(high) covariate <-  -covariate
    covcdf <- spatialcdf(covariate) 
    p <- seq(0, 1, length.out = nsteps)
    thresh <- quantile.ewcdf(covcdf, p)
    if(high) thresh <- -thresh
  } else {
    p      <- tp$p
    thresh <- tp$thresh
    if(is.null(p) || is.null(thresh))
      stop("tp should be a data frame with columns named 'p' and 'thresh'")
  }
  integrand <- X(thresh)
  est <- cumsum(integrand)
  est <- est/max(est)
  ## Basic data.frame
  df <- data.frame(p = p, null = p, thresh = thresh, est=est)
  ## Correct column name for estimate
  names(df)[4] <- if(method == "function") "fun" else method
  ## compute ROC 
  result <- as.roc.data.frame(df, method=method)
  return(result)
}


#' ROC methods for other classes

roc.spatialCDFframe <- function(X, ..., high=TRUE, plength=1024) {
  trap.extra.arguments(...)
  U <- X$values$U
  UU <- if(high) 1 - U else U
  ec <- ewcdf(UU, weights=X$weights)
  p <- seq(0, 1, length=plength)
  FZ <- X$values$FZ
  FZinverse <- quantilefun(FZ)
  thresh <- if(high) FZinverse(1-p) else FZinverse(p)
  #' set up data
  df <- data.frame(p=p, thresh=thresh, null=p, raw=ec(p))
  #' Convert to roc+fv object with nice labels etc.
  rslt <- as.roc.data.frame(df, method = "raw")
  return(rslt)
}

roc.cdftest <- function(X, ..., high=TRUE) {
  roc(attr(X, "frame"), high=high, ...)
}

roc.bermantest <- function(X, ..., high=TRUE) {
  roc(X[["fram"]], high=high, ...)
}


#'
#'   Compute estimate of ROC by kernel smoothing
#'   and compute confidence bands
#'
#'   Original code by Suman Rakshit 2017, edited by Adrian Baddeley 2018-2023
#' 

rocSmoothCalc <- function(ZX, ZU,
                          ...,
                          weightsX=NULL, weightsU=NULL,
                          high=TRUE,
                          kernel="gaussian", bw="nrd0", adjust=1, 
                          alpha=0.05,
                          nGrid=2048,
                          p=seq(0, 1, length.out=plength),
                          plength=1024,
                          doCI=TRUE,
                          degfreeU=length(ZU)) {
  ra <- range(ZX, ZU, na.rm=TRUE)
  if(!is.null(weightsX)) {
    check.nvector(weightsX, length(ZX), things="points of X")
    weightsX <- weightsX/sum(weightsX)
  }
  if(!is.null(weightsU)) {
    check.nvector(weightsU, length(ZU), things="points/pixels of U")
    weightsU <- weightsU/sum(weightsU)
  }
  if(!missing(p)) {
    check.nvector(p)
    stopifnot(min(p) >= 0)
    stopifnot(max(p) <= 1)
    if(prod(range(diff(p))) < 0) stop("p should be a monotone sequence")
    plength <- length(p)
  }
  #' bw and adjust can be vectors of length 2
  #' specifying smoothing for ROC [1] and densities [2]
  bw <- rep(bw, 2)
  adjust <- rep(adjust, 2)
  #' degfreeU is either a single number (including Inf) or a function
  if(is.null(degfreeU)) degfreeU <- length(ZU) 
  if(is.function(degfreeU)) degfreeU <- degfreeU(length(ZU))
  check.1.real(degfreeU)
  stopifnot(degfreeU > 0)
  #' transform observed values using empirical cdf of null
  uX <- if(high) ewcdf(-ZU, weightsU)(-ZX) else ewcdf(ZU, weightsU)(ZX)
  #' perform density estimation on this scale
  fuX <- density(uX, from=0, to=1, n=nGrid, bw=bw[1], adjust=adjust[1])
  #' integrate to obtain ROC estimate
  ppp <- fuX$x
  yyy <- fuX$y
  ROCfun <- approxfun(ppp, cumsum(yyy)/sum(yyy), yleft=0, yright=1, rule=1)
  ROCp <- ROCfun(p)
  #' estimate pdf's
  fX <- density(ZX, weights=weightsX,
                kernel=kernel, bw=bw[2], adjust=adjust[2],
                from=ra[1], to=ra[2], n=nGrid)
  sigma <- fX$bw
  fU <- density(ZU, weights=weightsU,
                kernel=kernel, bw=sigma,
                from=ra[1], to=ra[2], n=nGrid)
  zz <- fX$x
  fXz <- fX$y
  fUz <- fU$y
  fX <- approxfun(zz, fXz, yleft=0, yright=0)
  fU <- approxfun(zz, fUz, yleft=0, yright=0)
  #' integrate to obtain cdf's
  FXz <- cumsum(fXz)/sum(fXz)
  FUz <- cumsum(fUz)/sum(fUz)
  FX <- approxfun(zz, FXz, yleft=0, yright=1, rule=1)
  FU <- approxfun(zz, FUz, yleft=0, yright=1, rule=1)
  FUinverse <- approxfun(FUz, zz, rule=2)
  #' Compute threshold at each desired value of p
  zzp <- FUinverse(if(high) 1-p else p)
  #' assemble result
  df <- data.frame(p = p , thresh=zzp, smooth=ROCp)
  if(doCI) {
    #' Compute standard error
    se2A <- (1/length(ZX)) * ROCp * (1 - ROCp)
    se2B <- (1/degfreeU) * p    * (1-p)      * ((fX(zzp)/fU(zzp))^2)
    se <- sqrt(se2A + se2B)
    ## Compute the upper and lower pointwise confidence band
    crit <- stats::qnorm(1-alpha/2)
    hi <- pmin(1.0, ROCp + crit * se)
    lo <- pmax(0.0, ROCp - crit * se)
    ## Add to result
    df <- cbind(df, data.frame(se=se, lo = lo, hi = hi))
  }
  return(df)
}

