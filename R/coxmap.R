#'
#'  coxmap.R
#'
#' Implementation of
#'  T.F. Cox
#' A method for mapping the dense and sparse regions of a forest stand
#' Applied Statistics 28 (1979) 14-19
#'
#' $Revision: 1.3 $ $Date: 2026/03/11 05:10:11 $

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
  if(n > 30) {
    ## asymptotic N(0,1) approximation
    dn <- qnorm(alpha)
    sn <- qnorm(1 - alpha)
  } else {
    #' construct CDF of T
    an <- n * (n+1)/2
    kk <- 1:n
    cok <- ((-1)^(n-kk))* (kk^n)/(factorial(kk) * factorial(n-kk))
    if(all(is.finite(cok)) && max(abs(cok) < 1e12)) {
      ## numerically stable
      #' CDF(x) - alpha
      f <- function(x, a) { 1 - sum(cok * exp(-(an + x/b1n)/kk)) - a}
      #' find quantiles
      dn <- uniroot(f, c(b2n, 3), a=alpha)$root
      sn <- uniroot(f, c(b2n, 3), a=1-alpha)$root
    } else {
      ## fall back on N(0,1) approximation
      dn <- qnorm(alpha)
      sn <- qnorm(1 - alpha)
    }
  }
  #' compute T at every location
  for(k in 1:n) {
    Rk <- as.im(distfun(X, k=k), ...)
    Yk <- (lambda * pi) * Rk^2
    Sn <- if(k == 1) Yk else (Sn + Yk)
  }
  Tn <- b1n * Sn + b2n
  ## apply quantile cutoffs
  Z <- cut(Tn, breaks=c(-Inf, dn, sn, Inf),
           labels=c("clumped", "neither", "sparse"))
  return(Z)
}
