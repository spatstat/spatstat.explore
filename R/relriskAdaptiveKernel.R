#'
#'  relriskAdaptiveKernel.R
#'
#'  Estimation of relative risk using adaptive bandwidth kde
#' 
#'  $Revision: 1.5 $ $Date: 2026/09/19 09:19:20 $
#'

relriskAdaptiveKernel <- function(X, ...) UseMethod("relriskAdaptiveKernel")

relriskAdaptiveKernel.ppp <- function(X, bw=NULL, ...,
                                      adjust=1,
                                      at=c("pixels", "points"),
                                      weights = NULL,
                                      relative=FALSE, normalise=FALSE,
                                      casecontrol=TRUE, control=1, case) {
  stopifnot(is.ppp(X))
  if(is.NAobject(X)) return(NAobject("list"))
  if("se" %in% names(list(...)))
    stop("Standard errors are not supported in relriskAdaptiveKernel.ppp",
         call.=FALSE)
  ## resolve arguments
  control.given <- !missing(control)
  if(!control.given) control <- NULL
  case.given <- !missing(case) && !is.null(case)
  if(!case.given) case <- NULL
  at <- match.arg(at)
  ## evaluate numerical weights (multiple columns not allowed)
  weights <- pointweights(X, weights=weights, parent=parent.frame())
  ## determine bandwidth at each data point
  bw <- resolve.adaptive.bandwidths(unmark(X), bw=bw, ...,
                                    adjust=adjust, warn=TRUE)
  ## handle each column of marks
  stopifnot(is.marked(X))
  nc <- NCOL(marks(X))
  result <- if(nc == 1) {
              ## for efficiency
              rrakpppEngine(X, 
                          ..., 
                          bw          = bw,
                          at          = at,
                          weights     = weights,
                          relative    = relative,
                          normalise   = normalise,
                          casecontrol = casecontrol,
                          control     = control,
                          case        = case)
            } else {
              context <- paste("In column", 1:nc, "of marks:")
              MA <- list(...,
                         bw          = bw, 
                         at          = at,
                         weights     = weights,
                         relative    = relative,
                         normalise   = normalise,
                         casecontrol = casecontrol,
                         control     = control,
                         case        = case)
              mapply(rrakpppEngine,
                     X=unstack(X),
                     context = context,
                     MoreArgs=MA,
                     SIMPLIFY=FALSE)
            }
  return(result)
}


rrakpppEngine <- local({

  rrakpppEngine <- function(X, bw=NULL, ...,
                            at=c("pixels", "points"),
                            weights = NULL, 
                            relative=FALSE, normalise=FALSE,
                            casecontrol=TRUE, control=1, case=NULL,
                            context=NULL) {
    stopifnot(is.ppp(X))
    stopifnot(is.multitype(X))
    at <- match.arg(at)
    weighted <- !is.null(weights)
    control.given <- !is.null(control)
    case.given <- !is.null(case)
    if(!control.given) control <- 1
    ## 
    npts <- npoints(X)
    marx <- marks(X)
    imarks <- as.integer(marx)
    types <- levels(marx)
    ntypes <- length(types)
    if(ntypes == 1)
      stop(paste(context, "Data contains only one type of points"))
    ## 
    casecontrol <- casecontrol && (ntypes == 2)
    if((control.given || case.given) && !(casecontrol || relative)) {
      aa <- c("control", "case")[c(control.given, case.given)]
      nn <- length(aa)
      warning(paste(context,
                    ngettext(nn, "Argument", "Arguments"),
                    paste(sQuote(aa), collapse=" and "),
                    ngettext(nn, "was", "were"),
                    "ignored, because relative=FALSE and",
                    if(ntypes==2) "casecontrol=FALSE" else
                    "there are more than 2 types of points"))
    }
    ## initialise error report
    uhoh <- NULL
    ## prepare for analysis
    Y <- split(X)
    splitbw <- split(bw, marx)
    splitweights <- if(weighted) split(weights, marx) else rep(list(NULL), ntypes)
    ## threshold for 0/0
    tinythresh <- 8 * .Machine$double.eps
    ## .........................................
    ## compute intensity estimates for each type
    ## .........................................
    switch(at,
           pixels = {
             ## intensity estimates of each type
             Deach <- as.solist(mapply(densityAdaptiveKernel,
                                       X=Y,
                                       bw=splitbw,
                                       MoreArgs=list(...),
                                       SIMPLIFY=FALSE))
             ## compute intensity estimate for unmarked pattern
             Dall <- im.apply(Deach, sum, check=FALSE)
           },
           points = {
             ## intensity estimates of each type **at each data point**
             ## dummy variable matrix
             dumm <- matrix(0, npts, ntypes)
             dumm[cbind(seq_len(npts), imarks)] <- 1
             colnames(dumm) <- types
             dummweights <- if(weighted) dumm * weights else dumm
             Deach <- densityAdaptiveKernel(unmark(X), bw=bw, ...,
                                            weights=dummweights)
             ## compute intensity estimate for unmarked pattern
             Dall <- rowSums(Deach)
           })
    ## .........................................
    ## compute probabilities/risks
    ## .........................................
    if(ntypes == 2 && casecontrol) {
      if(control.given || !case.given) {
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:2)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol))
            stop(paste(context, "No points have mark =", sQuote(control)))
        } else
          stop(paste(context,
                     "Unrecognised format for argument", sQuote("control")))
        if(!case.given)
          icase <- 3 - icontrol
      }
      if(case.given) {
        stopifnot(length(case) == 1)
        if(is.numeric(case)) {
          icase <- case <- as.integer(case)
          stopifnot(case %in% 1:2)
        } else if(is.character(case)) {
          icase <- match(case, types)
          if(is.na(icase))
            stop(paste(context, "No points have mark =", sQuote(case)))
        } else stop(paste(context,
                          "Unrecognised format for argument", sQuote("case")))
        if(!control.given) 
          icontrol <- 3 - icase
      }
      ## compute ......
      switch(at,
             pixels = {
               ## compute probability of case
               Dcase <- Deach[[icase]]
               pcase <- Dcase/Dall
               ## correct small numerical errors
               pcase <- Clamp01(pcase)
               ## trap NaN values, and similar
               dodgy <- (Dall < tinythresh)
               nbg <- Badvalues(pcase) | Really(dodgy)
               if(any(nbg)) {
                 warning(paste(context,
                               "Numerical underflow detected:",
                               "sigma is probably too small"),
                         call.=FALSE)
                 uhoh <- unique(c(uhoh, "underflow"))
                 ## apply l'Hopital's rule:
                 ##     p(case) = 1{nearest neighbour is case}
                 distcase <- distmap(Y[[icase]], xy=pcase)
                 distcontrol <- distmap(Y[[icontrol]], xy=pcase)
                 closecase <- eval.im(as.integer(distcase < distcontrol))
                 pcase[nbg] <- closecase[nbg]
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- eval.im(ifelse(pcase < 1, pcase/(1-pcase), NA))
               }
             },
             points={
               ## compute probability of case
               pcase <- Deach[,icase]/Dall
               ## correct small numerical errors
               pcase <- Clamp01(pcase)
               ## trap NaN values
               dodgy <- (Dall < tinythresh)
               if(any(nbg <- Badvalues(pcase) | Really(dodgy))) {
                 warning(paste(context,
                               "Numerical underflow detected:",
                               "sigma is probably too small"),
                         call.=FALSE)
                 uhoh <- unique(c(uhoh, "underflow"))
                 ## apply l'Hopital's rule
                 nntype <- imarks[nnwhich(X)]
                 pcase[nbg] <- as.integer(nntype[nbg] == icase)
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- ifelse(pcase < 1, pcase/(1-pcase), NA)
               }
             })
    } else {
      ## several types
      if(relative) {
        ## need 'control' type
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:ntypes)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol))
            stop(paste(context, "No points have mark =", sQuote(control)))
        } else
          stop(paste(context,
                     "Unrecognised format for argument", sQuote("control")))
      }
      switch(at,
             pixels={
               #' Ops.imagelist not yet working
               probs <- imagelistOp(Deach, Dall, "/")
               ## correct small numerical errors
               probs <- as.solist(lapply(probs, Clamp01))
               ## trap NaN values
               nbg <- lapply(probs, Badvalues)
               nbg <- Reduce("|", nbg)
               dodgy <- (Dall < tinythresh)
               nbg <- nbg | Really(dodgy)
               if(any(nbg)) {
                 warning(paste(context,
                               "Numerical underflow detected:",
                               "sigma is probably too small"),
                         call.=FALSE)
                 uhoh <- unique(c(uhoh, "underflow"))
                 ## apply l'Hopital's rule
                 distX <- distmap(X, xy=Dall)
                 whichnn <- attr(distX, "index")
                 typenn <- eval.im(imarks[whichnn])
                 typennsub <- as.matrix(typenn)[nbg]
                 for(k in seq_along(probs)) 
                   probs[[k]][nbg] <- (typennsub == k)
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- as.solist(lapply(probs,
                                            Divideifpositive,
                                            d = probs[[icontrol]]))
               }
             },
             points = {
               probs <- Deach/Dall
               ## correct small numerical errors
               probs <- Clamp01(probs)
               ## trap NaN values
               dodgy <- (Dall < tinythresh)
               bad <- Badvalues(probs) 
               badrow <- matrowany(bad) | Really(dodgy)
               if(any(badrow)) {
                 warning(paste(context,
                               "Numerical underflow detected:",
                               "sigma is probably too small"),
                         call.=FALSE)
                 uhoh <- unique(c(uhoh, "underflow"))
                 ## apply l'Hopital's rule
                 typenn <- imarks[nnwhich(X)]
                 probs[badrow, ] <- (typenn == col(result))[badrow, ]
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- probs/probs[,icontrol]
               }
            })
    }
    attr(result, "bw") <- bw
    if(length(uhoh)) attr(result, "warnings") <- uhoh
    return(result)
  }

  Clamp01 <- function(x) {
    if(is.im(x)) return(eval.im(pmin(pmax(x, 0), 1)))
    return(pmin(pmax(x, 0), 1))
  }

  Badvalues <- function(x) {
    if(is.im(x)) x <- as.matrix(x)
    return(!(is.finite(x) | is.na(x)))
  }

  Really <- function(x) {
    if(is.im(x)) x <- as.matrix(x)
    x[is.na(x)] <- FALSE
    return(x)
  }
  
  Divideifpositive <- function(z, d) { eval.im(ifelse(d > 0, z/d, NA)) }
  
  rrakpppEngine
})


