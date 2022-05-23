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
#'  $Revision: 1.9 $ $Date: 2020/11/02 06:39:23 $

if(FULLTEST) {
local({

  hopskel.test(redwood, method="MonteCarlo", nsim=5)
  
  #' quadrat test - spatial methods
  a <- quadrat.test(redwood, 3)
  domain(a)
  shift(a, c(1,1))

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
#   $Revision: 1.35 $   $Date: 2022/04/14 00:49:39 $
#


