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
#
# tests/NAinCov.R
#
# Testing the response to the presence of NA's in covariates
#
# $Revision: 1.9 $ $Date: 2023/11/05 01:45:36 $

if(FULLTEST) {
local({
  #' quantile.ewcdf
  f <- ewcdf(runif(100), runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
  #' quantile.density
  f <- density(runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
})
}
#
#    tests/nnstat.R
#
# Check code that uses nndist/nnwhich
#
# nnorient()
# stienen()
#
#   $Revision: 1.1 $  $Date: 2020/12/04 03:45:44 $
#


local({
  if(FULLTEST) {
    #' test nnorient
    nnorient(cells, domain=erosion(Window(cells), 0.1))
    #' degenerate case
    X <- cells[nndist(cells) > bdist.points(cells)]
    f <- nnorient(X)
    #' nnclean
    A <- nnclean(shapley, k=17, edge.correct=TRUE)
    B <- nnclean(runifpoint3(300), 3)
    #' stienen set
    #' bug when disc radius is zero
    Y <- unmark(humberside)[40:100] # contains duplicated points
    stienen(Y)
    Z <- stienenSet(Y)
    #' other cases
    U <- stienen(cells[1])
    V <- stienenSet(cells, edge=FALSE)
  }
})


  
