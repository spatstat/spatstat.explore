#'
#'  coxmap.R
#'
#' Implementation of
#'  T.F. Cox
#' A method for mapping the dense and sparse regions of a forest stand
#' Applied Statistics 28 (1979) 14-19
#'
#' $Revision: 1.8 $ $Date: 2026/03/11 11:13:37 $

coxmap <- function(X, ...) {
  UseMethod("coxmap")
}

coxmap.ppp <- function(X, n, alpha=0.05, ...) {
  verifyclass(X, "ppp")
  check.1.integer(n)
  stopifnot(n > 1)
  lambda <- intensity(unmark(X))
  #' coefficients in statistic T
  b1n <- 1/sqrt(n * (n+1) * (2*n+1)/6)
  b2n <- -sqrt(3 * n * (n+1)/(2 * (2*n + 1)))
  #' find quantiles of T
  use.normal.approx <- (n >= 30)
  if(!use.normal.approx) {
    #' construct CDF of T
    an <- n * (n+1)/2
    kk <- 1:n
    cok <- ((-1)^(n-kk))* (kk^n)/(factorial(kk) * factorial(n-kk))
    #' check for numerical weirdness
    if(any(!is.finite(cok)) || (max(abs(cok)) > 1e12)) {
      use.normal.approx <- TRUE
    } else {
      ## numerically stable
      #' Find quantiles as roots of function f(x, alpha) = CDF(x) - alpha
      f <- function(x, a) { 1 - sum(cok * exp(-(an + x/b1n)/kk)) - a}
      dn <- try(uniroot(f, c(b2n, 3), a=alpha)$root, silent=TRUE)
      sn <- try(uniroot(f, c(b2n, 3), a=1-alpha)$root, silent=TRUE)
      if(inherits(dn, "try-error") || inherits(sn, "try-error")) {
        #' root-finding failed
        use.normal.approx <- TRUE
      }
    }
  }
  if(use.normal.approx) {
    ## asymptotic N(0,1) approximation
    dn <- qnorm(alpha)
    sn <- qnorm(1 - alpha)
  }
  #' ---------- calculate statistic T --------------------
  #' compute squared distances at every location
  R2 <- nnmap(X, k=1:n, what="dist", squared=TRUE, ...)
  #' sum them pixelwise
  Sn <- (lambda * pi) * im.apply(R2, sum)
  #' evaluate statistic
  Tn <- b1n * Sn + b2n
  #' --------- classify -------------------------------------
  ## apply quantile cutoffs
  Z <- cut(Tn, breaks=c(-Inf, dn, sn, Inf),
           labels=c("clumped", "neither", "sparse"))
  return(Z)
}
