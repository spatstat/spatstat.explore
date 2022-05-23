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
#'
#'     tests/threedee.R
#'
#'     Tests of 3D code 
#'
#'      $Revision: 1.8 $ $Date: 2020/05/02 01:32:58 $
#'

local({
  X <- runifpoint3(30)
  Y <- runifpoint3(20)
  if(FULLTEST) {
    A <- runifpoint3(10, nsim=2)
    Z <- ppsubset(X, 2:4)
  }
  ##
  if(ALWAYS) { # includes C code
    d <- pairdist(X, periodic=TRUE, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE, periodic=TRUE)
    #' 
    h <- has.close(X, 0.2)
    h <- has.close(X, 0.2, periodic=TRUE)
    h <- has.close(X, 0.2, Y=Y)
    h <- has.close(X, 0.2, Y=Y, periodic=TRUE)
    #' code blocks not otherwise reached
    rmax <- 0.6 * max(nndist(X))
    g <- G3est(X, rmax=rmax, correction="rs")
    g <- G3est(X, rmax=rmax, correction="km")
    g <- G3est(X, rmax=rmax, correction="Hanisch")
    g <- G3est(X, rmax=rmax, sphere="ideal")
    g <- G3est(X, rmax=rmax, sphere="digital")
    v <- sphere.volume()
    v <- digital.volume()
    #' older code
    co <- coords(X)
    xx <- co$x
    yy <- co$y
    zz <- co$z
    gg1 <- g3engine(xx, yy, zz, correction="Hanisch G3")
    gg2 <- g3engine(xx, yy, zz, correction="minus sampling")
    ff1 <- f3engine(xx, yy, zz, correction="no")
    ff2 <- f3engine(xx, yy, zz, correction="minus sampling")
  }
  ##
  if(ALWAYS) {
    #'class support
    X <- runifpoint3(10)
    print(X)
    print(X %mark% runif(10))
    print(X %mark% factor(letters[c(1:5,5:1)]))
    print(X %mark% data.frame(a=1:10, b=runif(10)))
    da <- as.Date(paste0("2020-01-0", c(1:5,5:1)))
    print(X %mark% da)
    print(X %mark% data.frame(a=1:10, b=da))
  }
})
