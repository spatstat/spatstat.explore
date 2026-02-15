#' resolve.lambda.R
#'
#' Common code to evaluate the intensity 'lambda' in Kinhom, pcfinhom
#'    (and multitype counterparts)
#'
#'    resolve.lambda
#'    resolve.lambdacross
#'    resolve.reciplambda
#'    validate.weights
#'    updateData (generic) (soon to be deprecated)
#'
#' $Revision: 1.26 $ $Date: 2026/02/15 09:17:24 $

resolve.lambda <- function(X, lambda=NULL, ...) {
  UseMethod("resolve.lambda")
}

resolve.lambda.ppp <- function(X, lambda=NULL, ...,
                               sigma=NULL, varcov=NULL,
                               leaveoneout=TRUE, update=TRUE,
                               check=TRUE) {
  dangerous <- "lambda"
  danger <- TRUE
  if(npoints(X) == 0) {
    danger <- FALSE
    lambda <- numeric(0)
  } else if(is.null(lambda)) {
    ## No intensity data provided
    ## Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=leaveoneout)
    lambda <- as.numeric(lambda)
    if(check) validate.weights(lambda, how="density estimation")
    danger <- FALSE
  } else if(is.im(lambda)) {
    lambda <- safelookup(lambda, X)
  } else if(is.function(lambda)) {
    lambda <- lambda(X$x, X$y)
  } else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
    check.nvector(lambda, npoints(X), vname="lambda")
  } else if(inherits(lambda, c("ppm", "kppm", "dppm", "slrm"))) {
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required when 'lambda' is a fitted model",
           call.=FALSE)
    ## model provides intensity
    model <- lambda
    if(update) {
      model <- update(model, X)
      danger <- FALSE
    }
    if(inherits(model, "slrm")) {
      #' predict.slrm has different syntax, 
      #' and does not support leave-one-out prediction
      lambda <- predict(model)[X]
    } else {
      lambda <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)

    }
  } else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image,",
                    "a fitted model, or a function"))
  if(check) validate.weights(lambda)
  return(list(lambda=lambda,
              danger=danger,
              dangerous=if(danger) dangerous else NULL))
}

resolve.reciplambda <- function(X, lambda=NULL, reciplambda=NULL, ...) {
  UseMethod("resolve.reciplambda")
}

resolve.reciplambda.ppp <- function(X, lambda=NULL, reciplambda=NULL,
                                    ...,
                                    sigma=NULL, varcov=NULL,
                                    leaveoneout=TRUE, update=TRUE,
                                    check=TRUE) {
  dangerous <- c("lambda", "reciplambda")
  danger <- TRUE
  if(npoints(X) == 0) {
    danger <- FALSE
    lambda <- reciplambda <- numeric(0)
  } else if(is.null(lambda) && is.null(reciplambda)) {
    ## No intensity data provided
    danger <- FALSE
    ## Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=leaveoneout)
    lambda <- as.numeric(lambda)
    if(check) validate.weights(lambda, how="density estimation")
    reciplambda <- 1/lambda
  } else if(!is.null(reciplambda)) {
    ## 1/lambda values provided
    lambda <- NULL
    if(is.im(reciplambda)) {
      reciplambda <- safelookup(reciplambda, X)
      if(check) validate.weights(reciplambda, recip=TRUE, how="image lookup")
    } else if(is.function(reciplambda)) {
      reciplambda <- reciplambda(X$x, X$y)
      if(check) validate.weights(reciplambda, recip=TRUE, how="function evaluation")
    } else if(is.numeric(reciplambda) && is.vector(as.numeric(reciplambda))) {
      check.nvector(reciplambda, npoints(X), vname="reciplambda")
      if(check) validate.weights(reciplambda, recip=TRUE)
    } else stop(paste(sQuote("reciplambda"),
                    "should be a vector, a pixel image, or a function"))
  } else {
    #' lambda values provided
    if(is.im(lambda)) {
      lambda <- safelookup(lambda, X)
      if(check) validate.weights(lambda, how="image lookup")
    } else if(is.function(lambda)) {
      lambda <- lambda(X$x, X$y)
      if(check) validate.weights(lambda, how="function evaluation")
    } else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
      check.nvector(lambda, npoints(X), vname="lambda")
      if(check) validate.weights(lambda)
    } else if(inherits(lambda, c("ppm", "kppm", "dppm", "slrm"))) {
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required when 'lambda' is a fitted model",
           call.=FALSE)
      ## model provides intensity
      model <- lambda
      if(update) {
        force(X)
        env.here <- sys.frame(sys.nframe())
        model <- update(model, X, envir=env.here)
        danger <- FALSE
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambda <- predict(model)[X]
      } else {
        lambda <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
      }
      if(check) validate.weights(lambda, how="model prediction")
    } else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image,",
                      "a fitted model, or a function"))
    ## evaluate reciprocal
    reciplambda <- 1/lambda
  }
  return(list(lambda=lambda, reciplambda=reciplambda,
              danger=danger,
              dangerous=if(danger) dangerous else NULL))
}

resolve.lambdacross <- function(X, I, J,
                                lambdaI=NULL, lambdaJ=NULL, ...) {
  UseMethod("resolve.lambdacross")
}

resolve.lambdacross.ppp <- function(X, I, J,
                                    lambdaI=NULL, lambdaJ=NULL,
                                    ...,
                                    lambdaX=NULL,
                                    sigma=NULL, varcov=NULL,
                                    leaveoneout=TRUE, update=TRUE,
                                    lambdaIJ=NULL,
                                    Iexplain="points satisfying condition I",
                                    Jexplain="points satisfying condition J",
                                    Ilevels=NULL,
                                    Jlevels=NULL) {
  dangerous <- c("lambdaI", "lambdaJ")
  dangerI <- dangerJ <- TRUE
  XI <- unmark(X[I])
  XJ <- unmark(X[J])
  nI <- npoints(XI)
  nJ <- npoints(XJ)
  nX <- npoints(X)

  lamIname <- short.deparse(substitute(lambdaI))
  lamJname <- short.deparse(substitute(lambdaJ))
  bothnames <- c(lamIname, lamJname)
  givenI <- !is.null(lambdaI)
  givenJ <- !is.null(lambdaJ)
  givenX <- !is.null(lambdaX)

  if(givenI != givenJ) {
    givenone <- bothnames[c(givenI, givenJ)]
    missedone <- setdiff(bothnames, givenone)
    stop(paste("If", givenone, "is given, then",
               missedone, "should also be given"),
         call.=FALSE)
  }
  if(givenX && givenI && givenJ)
    warning(paste(paste(sQuote(bothnames), collapse=" and "),
                  "were ignored because", sQuote("lambdaX"),
                  "was given"),
            call.=FALSE)

  if(givenX) {
    ## Single object must provide intensities for two processes X[I] and X[J]
    if(!is.multitype(X))
      stop("When lambdaX is given, X must be a multitype point pattern")
    if(is.null(Ilevels) || is.null(Jlevels))
      stop("When lambdaX is given, arguments Ilevels and Jlevels are required")
    Xlevels <- levels(marks(X))
    nM <- length(Xlevels)
    if(anyNA(match(c(Ilevels, Jlevels), Xlevels)))
      stop("Ilevels and Jlevels must be levels of the factor marks(X)")
    if(is.im(lambdaX)) {
      stop("A single image is not sufficient for lambdaX")
    } else if(is.vector(lambdaX)) {
      stop("A single vector is not sufficient for lambdaX")
    } else if(is.imlist(lambdaX) && length(lambdaX) == length(Xlevels)) {
      ## Look up intensity values
      lamXM <- sapply(lambdaX, "[", i=unmark(X), drop=FALSE)
      lambdaI <- rowSums(lamXM[I, match(Ilevels, Xlevels), drop=FALSE])
      lambdaJ <- rowSums(lamXM[J, match(Jlevels, Xlevels), drop=FALSE])
    } else if(is.function(lambdaX)) { 
      ## evaluate function f(x,y,m) for all locations, for each possible type
      lamXM <- matrix(, nX, nM) 
      for(m in 1:nM) 
        lamXM[,m] <- lambdaX(X$x, X$y, rep(Xlevels[m], nX))
      lambdaI <- rowSums(lamXM[I, match(Ilevels, Xlevels), drop=FALSE])
      lambdaJ <- rowSums(lamXM[J, match(Jlevels, Xlevels), drop=FALSE])
    } else if(is.numeric(lambdaX) && !is.null(dim(lambdaX))) {
      ## matrix of intensity values
      if(nrow(lambdaX) == nX && ncol(lambdaX) == nM) {
        lambdaI <- rowSums(lambdaX[I, match(Ilevels, Xlevels), drop=FALSE])
        lambdaJ <- rowSums(lambdaX[J, match(Jlevels, Xlevels), drop=FALSE])
      } else {
        stop(paste("Wrong dimensions for matrix lambdaX:",
                   paste(dim(lambdaX), collapse=" x "),
                   "should be",
                   paste(c(nX, nM), collapse=" x ")))
      }
    } else if(inherits(lambdaX, c("ppm", "kppm", "dppm"))) {
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required when 'lambdaX' is a fitted model",
           call.=FALSE)
      ## point process model provides intensity of each type of point
      model <- lambdaX
      if(!is.multitype(model))
        stop("Argument lambdaX should be a multitype point process model")
      if(update) {
        force(X)
        env.here <- sys.frame(sys.nframe())
        model <- update(model, X, envir=env.here)
        dangerI <- dangerJ <- FALSE
        dangerous <- "lambdaIJ"
      }
      ## predict model at observed locations for each possible mark
      lamXM <- matrix(, nX, nM) 
      for(m in 1:nM) 
        lamXM[,m] <- predict(model, locations = X %mark% Xlevels[m])
      ## overwrite values for observed locations and *observed* marks
      ## using leave-one-out predictions
      if(update) {
        fitlam <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        if(length(fitlam) != nX)
          stop(paste("Internal error: trying to overwrite",
                     nX, "predicted values with",
                     length(fitlam), "fitted values"),
               call.=FALSE)
        lamXM[cbind(1:nX, as.integer(marks(X)))] <- fitlam
      }
      ## now compute marginal intensities
      lambdaI <- rowSums(lamXM[I, match(Ilevels, Xlevels), drop=FALSE])
      lambdaJ <- rowSums(lamXM[J, match(Jlevels, Xlevels), drop=FALSE])
    } else if(inherits(lambdaX, "slrm")) {
      ## does not (yet) support multitype processes
      stop("lambdaX cannot be a model of class 'slrm'")
    } else stop(paste("Argument lambdaX is not understood:",
                      "it should be a numeric matrix,",
                      "a list of images, a function(x,y,marks)",
                      "or a fitted multitype point process model (ppm, kppm or dppm)"))
  } else {
    ## lambdaI, lambdaJ expected
    if(!givenI) {
      ## estimate intensity
      dangerI <- FALSE
      dangerous <- setdiff(dangerous, "lambdaI")
      lambdaI <- density(XI, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaI)) {
      ## look up intensity values
      lambdaI <- safelookup(lambdaI, XI)
    } else if(is.function(lambdaI)) {
      ## evaluate function at locations of X[I]
      lambdaI <- lambdaI(XI$x, XI$y)
    } else if(is.numeric(lambdaI) && is.vector(as.numeric(lambdaI))) {
      ## validate intensity vector
      check.nvector(lambdaI, nI, things=Iexplain, vname="lambdaI")
    } else if(inherits(lambdaI, c("ppm", "kppm", "dppm", "slrm"))) {
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required when 'lambdaI' is a fitted model",
           call.=FALSE)
      ## point process model provides intensity
      model <- lambdaI
      if(update) {
        env.here <- sys.frame(sys.nframe())
        if(is.multitype(model)) {
          ## multitype model provides intensity of each type of point
          ## assumed fitted to X
          force(X)
          model <- update(model, X, envir=env.here)
        } else {
          ## unmarked model provides intensity of subset I
          ## assumed fitted to X[I]
          model <- update(model, XI, envir=env.here)
        }
        dangerI <- FALSE
        dangerous <- setdiff(dangerous, "lambdaI")
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambdaI <- predict(model)[XI]
      } else if(!is.multitype(model)) {
        ## unmarked model assumed to have been fitted to X[I]
        if(update) {
          lambdaI <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        } else {
          lambdaI <- predict(model, locations=XI)
        }
      } else {
        ## multitype model assumed fitted to X
        if(is.null(Ilevels))
          stop(paste("When lambdaI is a multitype point process,",
                     "argument Ilevels is required"))
        if(!is.multitype(X))
          stop(paste("When lambdaI is a multitype point process,",
                     "X should be multitype"))
        Xlevels <- levels(marks(X))
        nM <- length(Xlevels)
        if(anyNA(match(Ilevels, Xlevels)))
          stop("Ilevels must contain levels of the factor marks(X)")
        lamXM <- matrix(, nX, nM) 
        for(m in 1:nM) 
          lamXM[,m] <- predict(model, locations = X %mark% Xlevels[m])
        ## overwrite values for observed locations and *observed* marks
        ## leave-one-out predictions
        if(update) {
          fitlam <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
          if(length(fitlam) != nX)
            stop(paste("Internal error: trying to overwrite",
                       nX, "predicted values with",
                       length(fitlam), "fitted values"),
                 call.=FALSE)
          lamXM[cbind(1:nX, as.integer(marks(X)))] <- fitlam
        }
        ## sum to get intensity of subset I
        lambdaI <- rowSums(lamXM[I, match(Ilevels, Xlevels), drop=FALSE])
      }
    } else stop(paste(sQuote("lambdaI"), "should be a vector or an image"))
    
    if(!givenJ) {
      ## estimate intensity
      dangerJ <- FALSE
      dangerous <- setdiff(dangerous, "lambdaJ")
      lambdaJ <- density(XJ, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaJ)) {
      ## look up intensity values
      lambdaJ <- safelookup(lambdaJ, XJ)
    } else if(is.function(lambdaJ)) {
      ## evaluate function at locations
      lambdaJ <- lambdaJ(XJ$x, XJ$y)
    } else if(is.numeric(lambdaJ) && is.vector(as.numeric(lambdaJ))) {
      ## validate intensity vector
      check.nvector(lambdaJ, nJ, things=Jexplain, vname="lambdaJ")
    } else if(inherits(lambdaJ, c("ppm", "kppm", "dppm", "slrm"))) {
    if(!requireNamespace("spatstat.model")) 
      stop("The package spatstat.model is required when 'lambdaJ' is a fitted model",
           call.=FALSE)
      ## point process model provides intensity
      model <- lambdaJ
      if(update) {
        env.here <- sys.frame(sys.nframe())
        force(X)
        if(is.multitype(model)) {
          ## multitype model provides intensity of each type of point
          ## assumed fitted to X
          model <- update(model, X, envir=env.here)
        } else {
          ## unmarked model provides intensity of subset I
          ## assumed fitted to X[J]
          model <- update(model, XJ, envir=env.here)
        }
        dangerJ <- FALSE
        dangerous <- setdiff(dangerous, "lambdaJ")
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambdaJ <- predict(model)[XJ]
      } else if(!is.multitype(model)) {
        ## unmarked model assumed to have been fitted to X[J]
        if(update) {
          lambdaJ <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        } else {
          lambdaJ <- predict(model, locations=XJ)
        }
      } else {
        ## multitype model assumed fitted to X
        if(is.null(Jlevels))
          stop(paste("When lambdaJ is a multitype point process,",
                     "argument Jlevels is required"))
        if(!is.multitype(X))
          stop(paste("When lambdaJ is a multitype point process,",
                     "X should be multitype"))
        Xlevels <- levels(marks(X))
        nM <- length(Xlevels)
        if(anyNA(match(Jlevels, Xlevels)))
          stop("Ilevels must contain levels of the factor marks(X)")
        lamXM <- matrix(, nX, nM) 
        for(m in 1:nM) 
          lamXM[,m] <- predict(model, locations = X %mark% Xlevels[m])
        ## overwrite values for observed locations and *observed* marks
        ## leave-one-out predictions
        if(update) {
          fitlam <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
          if(length(fitlam) != nX)
            stop(paste("Internal error: trying to overwrite",
                       nX, "predicted values with",
                       length(fitlam), "fitted values"),
                 call.=FALSE)
          lamXM[cbind(1:nX, as.integer(marks(X)))] <- fitlam
        }
        ## sum to get intensity of subset J
        lambdaJ <- rowSums(lamXM[J, match(Jlevels, Xlevels), drop=FALSE])
      }
    } else 
      stop(paste(sQuote("lambdaJ"), "should be a vector or an image"))
  }
  
  ## Product of weights for each pair
  if(!is.null(lambdaIJ)) {
    dangerIJ <- TRUE
    dangerous <- union(dangerous, "lambdaIJ")
    if(!is.matrix(lambdaIJ))
      stop("lambdaIJ should be a matrix")
    if(nrow(lambdaIJ) != nI)
      stop(paste("nrow(lambdaIJ) should equal the number of", Iexplain))
    if(ncol(lambdaIJ) != nJ)
      stop(paste("ncol(lambdaIJ) should equal the number of", Jexplain))
  } else {
    dangerIJ <- FALSE
  }
    
  danger <- dangerI || dangerJ || dangerIJ
    
  return(list(lambdaI = lambdaI, lambdaJ = lambdaJ, lambdaIJ=lambdaIJ,
                danger = danger, dangerous = dangerous))
}


validate.weights <- function(x, recip=FALSE, how = NULL,
                             allowzero = recip,
                             allowinf  = !recip) {
  ra <- range(x)
  offence <-
    if(!allowinf && !all(is.finite(ra)))  "infinite" else
    if(ra[1] < 0)                         "negative" else
    if(!allowzero && ra[1] == 0)          "zero" else NULL
  if(!is.null(offence)) {
    xname <- deparse(substitute(x))
    offenders <- paste(offence, "values of", sQuote(xname))
    if(is.null(how))
      stop(paste(offenders, "are not allowed"), call.=FALSE)
    stop(paste(how, "yielded", offenders), call.=FALSE)
  }
  return(TRUE)
}

## The following internal functions will soon be removed.
## They will be replaced by the idiom
##      model <- update(model, X, ...)

updateData <- function(model, X, ...) {
  UseMethod("updateData")
}

updateData.default <- function(model, X, ..., warn=TRUE) {
  ## for some bizarre reason, method dispatch often fails for this function
  ## so we do it by hand as a backup
  if(inherits(model, c("ppm", "kppm", "dppm", "slrm"))) {
    if (requireNamespace("spatstat.model")) {
      model <- update(model, X)
    } else if (warn) {
      warning("Model was not updated; this requires package spatstat.model", 
              call. = FALSE)
    }
  } else if(inherits(model, "lppm")) {
      if (requireNamespace("spatstat.linnet")) {
        model <- update(model, X)
      } else if (warn) {
        warning("Model was not updated; this requires package spatstat.linnet", 
                call. = FALSE)
      }
  } else if (warn) {
    warning("Unrecognised kind of 'model'; no update performed", 
            call. = FALSE)
  }
  return(model)
}

