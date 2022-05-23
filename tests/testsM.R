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
##    tests/markcor.R
##
##   Tests of mark correlation code (etc)
##
## $Revision: 1.7 $ $Date: 2020/11/25 01:23:32 $

local({
  if(ALWAYS) {
    ## check.testfun checks equality of functions
    ##  and is liable to break if the behaviour of all.equal is changed
    fe <- function(m1, m2) {m1 == m2}
    fm <- function(m1, m2) {m1 * m2}
    fs <- function(m1, m2) {sqrt(m1)}
    if(check.testfun(fe, X=amacrine)$ftype != "equ")
      warning("check.testfun fails to recognise mark equality function")
    if(check.testfun(fm, X=longleaf)$ftype != "mul")
      warning("check.testfun fails to recognise mark product function")
    check.testfun(fs, X=longleaf)
    check.testfun("mul")
    check.testfun("equ")
  }

  if(FULLTEST) {
    ## test all is well in Kmark -> Kinhom 
    MA <- Kmark(amacrine,function(m1,m2){m1==m2})
    set.seed(42)
    AR <- rlabel(amacrine)
    MR <- Kmark(AR,function(m1,m2){m1==m2})
    if(isTRUE(all.equal(MA,MR)))
      stop("Kmark unexpectedly ignores marks")

    ## cover code blocks in markcorr()
    X <- runifpoint(100) %mark% runif(100)
    Y <- X %mark% data.frame(u=runif(100), v=runif(100))
    ww <- runif(100)
    fone <- function(x) { x/2 }
    ffff <- function(x,y) { fone(x) * fone(y) }
    aa <- markcorr(Y)
    bb <- markcorr(Y, ffff, weights=ww, normalise=TRUE)
    bb <- markcorr(Y, ffff, weights=ww, normalise=FALSE)
    bb <- markcorr(Y, f1=fone, weights=ww, normalise=TRUE)
    bb <- markcorr(Y, f1=fone, weights=ww, normalise=FALSE)

    ## markcrosscorr
    a <- markcrosscorr(betacells, normalise=FALSE)
    if(require(sm)) {
      b <- markcrosscorr(betacells, method="sm")
    }

    ## Vmark with normalisation
    v <- Vmark(spruces, normalise=TRUE)
    v <- Vmark(finpines, normalise=TRUE)
  }
})
#' tests/mctests.R
#' Monte Carlo tests
#'        (mad.test, dclf.test, envelopeTest, hasenvelope)
#' $Revision: 1.5 $ $Date: 2022/05/23 04:09:49 $

local({
  if(FULLTEST) {
    envelopeTest(cells, Lest, exponent=1, nsim=9, savepatterns=TRUE)
    (a3 <- envelopeTest(cells, Lest, exponent=3, nsim=9, savepatterns=TRUE))
    
    envelopeTest(a3, Lest, exponent=3, nsim=9, alternative="less")
    
    envelopeTest(redwood, Lest, exponent=1, nsim=19,
                 rinterval=c(0, 0.1), alternative="greater", clamp=TRUE)
    envelopeTest(redwood, pcf, exponent=Inf, nsim=19,
                 rinterval=c(0, 0.1), alternative="greater", clamp=TRUE)
  }
})


