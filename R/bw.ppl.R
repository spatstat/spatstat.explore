#'
#'   bw.ppl.R
#'
#'   Likelihood cross-validation for kernel smoother of point pattern
#'
#'   $Revision: 1.17 $ $Date: 2024/09/06 07:09:13 $
#'

bw.ppl <- function(X, ..., srange=NULL, ns=16, sigma=NULL, varcov1=NULL, 
                   weights=NULL, shortcut=TRUE, warn=TRUE) {
  stopifnot(is.ppp(X))
  if(!is.null(varcov1))
    check.nmatrix(varcov1, 2, things="spatial dimensions", mname="varcov1")
  if(!is.null(sigma)) {
    stopifnot(is.numeric(sigma) && is.vector(sigma))
    ns <- length(sigma)
  } else {
    if(!is.null(srange)) {
      check.range(srange)
    } else {
      ## default rule based on point pattern spacing and window size
      nnd <- nndist(X)
      srange <- c(min(nnd[nnd > 0]), diameter(as.owin(X))/2)
      if(!is.null(varcov1)) {
        dref <- det(varcov1)^(1/4)
        srange <- srange/dref
      }
    }
    sigma <- geomseq(from=srange[1L], to=srange[2L], length.out=ns)
  }
  cv <- numeric(ns)
  if(shortcut) {
    for(i in 1:ns) {
      if(is.null(varcov1)) {
        si <- sigma[i]
        vi <- NULL
      } else {
        si <- NULL
        vi <- (sigma[i]^2) * varcov1
      }
      lamx <- density(X, sigma=si, varcov=vi,
                      at="points", leaveoneout=TRUE,
                      weights=weights, ...)
      cv[i] <- sum(log(lamx))
    }
  } else {
    IntLam <- numeric(ns)
    for(i in 1:ns) {
      if(is.null(varcov1)) {
        si <- sigma[i]
        vi <- NULL
      } else {
        si <- NULL
        vi <- (sigma[i]^2) * varcov1
      }
      lamx <- density(X, sigma=si, varcov=vi,
                      at="points", leaveoneout=TRUE,
                      weights=weights, ...)
      lam <- density(X, sigma=si, varcov=vi,
                     weights=weights, ...)
      mu <- integral.im(lam)
      cv[i] <- sum(log(lamx)) - mu
      IntLam[i] <- mu
    }
  }
  result <- bw.optim(cv, sigma, iopt=which.max(cv),
                     optimum="max",
                     creator="bw.ppl",
                     criterion="Likelihood Cross-Validation",
                     warnextreme=warn,
                     hargnames="srange",
                     unitname=if(is.null(varcov1)) unitname(X) else NULL,
                     template=varcov1, exponent=2)
  if(!shortcut) 
    attr(result, "info") <- list(IntegralLambda=IntLam)
  return(result)
}


