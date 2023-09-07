#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.explore
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.explore)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
##
##    tests/gcc323.R
##
##    $Revision: 1.3 $  $Date: 2020/04/28 12:58:26 $
##
if(ALWAYS) { # depends on hardware
local({
  # critical R values that provoke GCC bug #323
  a <- marktable(lansing, R=0.25)
  a <- marktable(lansing, R=0.21)
  a <- marktable(lansing, R=0.20)
  a <- marktable(lansing, R=0.10)
})
}
#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.10 $ $Date: 2023/07/17 07:30:48 $

if(FULLTEST) {
local({

  hopskel.test(redwood, method="MonteCarlo", nsim=5)
  
  #' quadrat test - spatial methods
  a <- quadrat.test(redwood, 3)
  domain(a)
  shift(a, c(1,1))

  #' quadrat test - correctness of mapping from table to quadrats
  Q2 <- quadratcount(humberside, 2, 3)
  T2 <- suppressWarnings(quadrat.test(Q2))
  R2 <- cbind(as.numeric(t(Q2)), round(10 * residuals(T2)))
  R2correct <- cbind(c(2, 20, 13, 11, 34, 123),
                     c(-46, -12, -62, -41, 50, 134))
  if(!all(R2 == R2correct))
    stop("Incorrect count-residual map for quadrat.test(2,3)")

  Q5 <- quadratcount(humberside, 5, 3)
  T5 <- suppressWarnings(quadrat.test(Q5))
  R5 <- cbind(as.numeric(t(Q5)), round(10 * residuals(T5)))
  R5correct <- cbind(
    c(  0,   0,   3, 19,   3,   2,  14,   5,  0,   2, 117, 35,   3),
    c(-19, -33, -42, 16, -37, -49, -28, -35, -5, -21, 295, 40, -32))
  if(!all(R5 == R5correct))
    stop("Incorrect count-residual map for quadrat.test(5,3)")
    
  #' cases of studpermu.test
  #' X is a hyperframe
  b <- studpermu.test(pyramidal, nperm=9)
  b <- studpermu.test(pyramidal, nperm=9, use.Tbar=TRUE)
  #' X is a list of lists of ppp
  ZZ <- split(pyramidal$Neurons, pyramidal$group)
  bb <- studpermu.test(ZZ, nperm=9)

  #' Issue #115
  X <- runifpoint(50, nsim = 3)
  Y <- runifpoint(3000, nsim = 3)
  h <- hyperframe(ppp = c(X, Y), group = rep(1:2, 3))
  studpermu.test(h, ppp ~ group)

  #' scan test
  Z <- scanmeasure(cells, 0.1, method="fft")
  rr <- c(0.05, 1)
  scan.test(amacrine, rr, nsim=5,
            method="binomial", alternative="less")
})
}
#
#  tests/imageops.R
#
#   $Revision: 1.43 $   $Date: 2023/08/29 01:03:59 $
#


if(FULLTEST) {
  local({
    #' case of "[.im" and "[<-.im" where index is an ssf
    d <- distmap(cells, dimyx=32)
    Empty <- cells[FALSE]
    EmptyFun <- ssf(Empty, numeric(0))
    ff <- d[EmptyFun]
    d[EmptyFun] <- 42

    #' Smooth.im -> blur.im with sigma=NULL
    Z <- as.im(function(x,y) { x - y }, letterR, dimyx=32)
    ZS <- Smooth(Z)

    #' deprecated -> im.apply(DA, which.max)
    Z <- which.max.im(bei.extra) 

    #' rotmean
    U <- rotmean(Z, origin="midpoint", result="im", padzero=FALSE)
    
    #' cases of distcdf
    distcdf(cells[1:5])
    distcdf(W=cells[1:5], dW=1:5)
    distcdf(W=Window(cells), V=cells[1:5])
    distcdf(W=Window(cells), V=cells[1:5], dV=1:5)
  })
}

