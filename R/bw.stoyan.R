#' bw.stoyan.R
#'
#' Stoyan's and Fiksel's rules for bandwidth for pcf estimation
#' with optional extension proposed by Baddeley, Davies and Hazelton
#'
#' D. Stoyan and H. Stoyan (1995)
#' Fractals, random shapes and point fields
#' Wiley, Chichester
#' 
#' T. Fiksel
#' Edge-corrected density estimators for point processes
#' Statistics 19 (1988) #1, 77-86
#' 
#' A. Baddeley, T.M.Davies and M.L.Hazelton
#' An improved estimator of the pair correlation function of a spatial point process.
#' Biometrika 111 (2025) #2, asaf021
#' 
#' Copyright (c) 2010-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#' $Revision: 1.3 $ $Date: 2026/04/30 06:50:25 $

bw.stoyan <- function(X, co = 0.15, extrapolate=FALSE, ...) 
{
  stopifnot(is.ppp(X))
  n <- max(1, npoints(X))
  W <- Window(X)
  a <- area(W)
  ## Stoyan rule for halfwidth of Epanechnikov kernel
  ## Stoyan & Stoyan 1995, eq (15.16), page 285
  hepa <- co/sqrt(n/a)
  ## convert to bandwidth (standard deviation)
  bw <- hepa/sqrt(5)
  ## extrapolate? (Baddeley Davies and Hazelton 2025)
  if(extrapolate)
    bw <- bw * ((100/n)^(1/5))
  return(bw)
}

bw.fiksel <- function(X, co = 0.1, extrapolate=FALSE, ...) 
{
  stopifnot(is.ppp(X))
  n <- max(1, npoints(X))
  W <- Window(X)
  a <- area(W)
  ## Fiksel rule for bandwidth (standard deviation) of Epanechnikov kernel
  bw <- co/sqrt(n/a)
  ## extrapolate? (Baddeley Davies and Hazelton 2025)
  if(extrapolate)
    bw <- bw * ((100/n)^(1/5))
  return(bw)
}

