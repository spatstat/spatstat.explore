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
## tests/rhohat.R
##
## Test all combinations of options for rhohatCalc
##
## $Revision: 1.6 $ $Date: 2022/05/22 08:03:48 $

local({
  if(FULLTEST) {
    X <-  rpoispp(function(x,y){exp(3+3*x)})
    Z <- as.im(function(x,y) { x }, Window(X))
    f <- funxy(function(x,y) { y + 1 }, Window(X))
    
    ## rhohat.ppp
    ## done in example(rhohat):
    ## rhoA <- rhohat(X, "x")
    ## rhoB <- rhohat(X, "x", method="reweight")
    ## rhoC <- rhohat(X, "x", method="transform")
    ## alternative smoother (if package locfit available)
    rhoA <- rhohat(X, "x", smoother="local")
    rhoB <- rhohat(X, "x", smoother="local", method="reweight")
    rhoC <- rhohat(X, "x", smoother="local", method="transform")

    #' code blocks
    rhoD <- rhohat(X, "y", positiveCI=TRUE)
    rhoE <- rhohat(X, Z,   positiveCI=TRUE)
    #' weights 
    rhoF <- rhohat(X, Z,   weights=f(X))
    rhoG <- rhohat(X, Z,   weights=f)
    rhoH <- rhohat(X, Z,   weights=as.im(f))
    
    lam <- as.im(function(x,y) {exp(3+2*x)}, W=Window(Z))

    ## Baseline
    rhoAb <- rhohat(X, "x", baseline=lam)
    rhoBb <- rhohat(X, "x", method="reweight", baseline=lam)
    rhoCb <- rhohat(X, "x", method="transform", baseline=lam)

    ## Horvitz-Thompson
    rhoAH <- rhohat(X, "x", horvitz=TRUE) 
    rhoBH <- rhohat(X, "x", method="reweight", horvitz=TRUE)
    rhoCH <- rhohat(X, "x", method="transform", horvitz=TRUE)

    ## class support
    plot(rhoA)
    plot(rhoA, rho ~ x, shade=NULL)
    plot(rhoA, log(rho) ~ x, shade=NULL)
    plot(rhoA, log(.) ~ x)

    ## rho2hat
    r2xy <- rho2hat(X, "x", "y")
    r2xyw <- rho2hat(X, "x", "y", method="reweight")
    print(r2xyw)
    plot(r2xy, do.points=TRUE)
    xcoord <- function(x,y) x
    ycoord <- function(x,y) y
    xim <- as.im(xcoord, W=Window(X))
    r2fi <- rho2hat(X, ycoord, xim)
    r2if <- rho2hat(X, xim, ycoord)
  }
})
