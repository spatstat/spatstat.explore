#'
#'  coxmap.R
#'
#' Implementation of
#'  T.F. Cox
#' A method for mapping the dense and sparse regions of a forest stand
#' Applied Statistics 28 (1979) 14-19
#'
#' $Revision: 1.1 $ $Date: 2026/01/21 06:26:39 $

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
  #' construct CDF of T
  kk <- 1:n
  cok <- ((-1)^(n-kk))* (kk^n)/(factorial(kk) * factorial(n-kk))
  an <- n * (n+1)/2
  #' CDF(x) - alpha
  f <- function(x, a) { 1 - sum(cok * exp(-(an + x/b1n)/kk)) - a}
  #' find quantiles
  dn <- uniroot(f, c(b2n, 3), a=alpha)$root
  sn <- uniroot(f, c(b2n, 3), a=1-alpha)$root
  #' compute T at every location
  for(k in 1:n) {
    Rk <- as.im(distfun(X, k=k), ...)
    Yk <- Rk^2
    Tn <- if(k == 1) Yk else (Tn + Yk)
  }
  Tn <- lambda * pi * Tn
  Tn <- b1n * Tn + b2n
  ## apply quantile cutoffs
  Z <- cut(Tn, breaks=c(-Inf, dn, sn, Inf),
           labels=c("clumped", "neither", "sparse"))
  return(Z)
}
