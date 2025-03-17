#'
#'   bw.bdh.R
#'
#'   'Rule of toe' for bandwidth
#'
#' Copyright (c) 2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#'  $Revision: 1.1 $ $Date: 2025/03/15 10:01:48 $

bw.bdh <- function(X, lambda=NULL, ..., base=bw.stoyan, k=2) {
  check.1.integer(k)
  stopifnot(k > 0)
  #' evaluate base bandwidth using rule 
  if(is.character(base)) base <- get(base, mode="function")
  if(is.function(base)) {
    bw <- do.call.matched(base, list(X, ...), matchfirst=TRUE)
    bw <- as.numeric(bw)
    check.bandwidth(bw, "value returned by 'base'")
  } else if(check.bandwidth(base, fatal=FALSE)) {
    bw <- base
  } else {
    stop("'base' should be a function or a single numeric value > 0",
         call.=FALSE)
  }
  #' evaluate intensity at points of X
  v <- resolve.reciplambda(X, lambda, ...)
  lam <- v$lambda
  mal <- v$reciplambda
  #' avoid divide-by-zero etc
  ok <- (lam > 0) & is.finite(mal)
  if(!all(ok)) {
    lam <- lam[ok]
    mal <- mal[ok]
  }
  #' go
  a <- mean(lam[ok]) * ((mean(mal[ok]^k))^(1/k))
  bw <- a * bw
  attr(bw, "adjust") <- a
  return(bw)
}

