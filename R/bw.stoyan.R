#' bw.stoyan.R
#'
#' Stoyan's rule for bandwidth for pcf estimation
#' with optional extension proposed by Baddeley, Davies and Hazelton
#' 
#' Copyright (c) 2010-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#' $Revision: 1.2 $ $Date: 2025/03/15 09:44:28 $

bw.stoyan <- function(X, co = 0.15, extrapolate=FALSE, ...) 
{
  stopifnot(is.ppp(X))
  n <- max(1, npoints(X))
  W <- Window(X)
  a <- area(W)
  bw <- co/sqrt(5 * n/a)
  if(extrapolate)
    bw <- bw * ((100/n)^(1/5))
  return(bw)
}

