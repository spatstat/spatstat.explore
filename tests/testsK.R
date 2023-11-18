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
#'  tests/kernelstuff.R
#'
#'  $Revision: 1.2 $  $Date: 2023/11/05 01:49:45 $

local({
  if(FULLTEST) {
    #' test all cases in kernels.R
    kernames <- c("gaussian", "rectangular", "triangular",
                  "epanechnikov", "biweight", "cosine", "optcosine")
    X <- rnorm(20)
    U <- runif(20)
    for(ker in kernames) {
      dX <- dkernel(X, ker)
      fX <- pkernel(X, ker)
      qU <- qkernel(U, ker)
      m0 <- kernel.moment(0, 0, ker)
      m1 <- kernel.moment(1, 0, ker)
      m2 <- kernel.moment(2, 0, ker)
      m3 <- kernel.moment(3, 0, ker)
    }
  }
})

#'
#'   tests/Kfuns.R
#'
#'   Various K and L functions and pcf
#'
#'   $Revision: 1.43 $  $Date: 2022/06/17 01:47:08 $
#'
#'   Assumes 'EveryStart.R' was run

if(FULLTEST) {
  Cells <- cells
  Amacrine <- amacrine
  Redwood <- redwood
} else {
  ## reduce numbers of data + dummy points
  spatstat.options(npixel=32, ndummy.min=16)
  Cells <- cells[c(FALSE, TRUE)]
  Amacrine <- amacrine[c(FALSE, TRUE)]
  Redwood <- redwood[c(FALSE, TRUE)]
}


myfun <- function(x,y){(x+1) * y } # must be outside

local({
  if(FULLTEST) {
    #' supporting code
    rmax.rule("Kscaled", owin(), 42)
    implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                      "polygonal", TRUE)
    implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                      "mask", TRUE)
    implemented.for.K(c("border", "isotropic"), "mask", TRUE)
    implemented.for.K(c("border", "isotropic"), "mask", FALSE)
    #' shortcuts
    D <- density(Cells)
    K <- Kborder.engine(Cells, rmax=0.4, weights=D, ratio=TRUE)
    K <- Knone.engine(Cells, rmax=0.4, weights=D, ratio=TRUE)
    allcor <- c("none", "border", "bord.modif","isotropic", "translate")
    K <- Krect.engine(Cells, rmax=0.4, ratio=TRUE, correction=allcor)
    K <- Krect.engine(Cells, rmax=0.4, ratio=TRUE, correction=allcor,
                      weights=D)
    K <- Krect.engine(Cells, rmax=0.4, ratio=TRUE, correction=allcor,
                      use.integers=FALSE)
    #' Kest special code blocks
    K <- Kest(Cells, var.approx=TRUE, ratio=FALSE)
    Z <- distmap(Cells) + 1
    Kb <- Kest(Cells, correction=c("border","bord.modif"),
               weights=Z, ratio=TRUE)
    Kn <- Kest(Cells, correction="none",
               weights=Z, ratio=TRUE)
    Knb <- Kest(Cells, correction=c("border","bord.modif","none"),
                weights=Z, ratio=TRUE)
  }
  if(ALWAYS) {
    bigint <- 50000 # This is only "big" on a 32-bit system where
                    # sqrt(.Machine$integer.max) = 46340.9
    X <- runifpoint(bigint)
    Z <- as.im(1/bigint, owin())
    Kb <- Kest(X, correction=c("border","bord.modif"),
               rmax=0.02, weights=Z, ratio=TRUE)
  }
  if(FULLTEST) {
    Kn <- Kest(X, correction="none",
               rmax=0.02, weights=Z, ratio=TRUE)
    Knb <- Kest(X, correction=c("border","bord.modif","none"),
                rmax=0.02, weights=Z, ratio=TRUE)
    #' pcf.ppp special code blocks
    pr  <- pcf(Cells, ratio=TRUE, var.approx=TRUE)
    pc  <- pcf(Cells, domain=square(0.5))
    pcr <- pcf(Cells, domain=square(0.5), ratio=TRUE)
    pw <- pcf(Redwood, correction="none")
    pwr <- pcf(Redwood, correction="none", ratio=TRUE)
    pv <- pcf(Redwood, kernel="rectangular")
    p1 <- pcf(Redwood[1])
    #' pcf.fv
    K <- Kest(Redwood)
    g <- pcf(K, method="a")
    g <- pcf(K, method="c")
    g <- pcf(K, method="d")
    #' Kinhom code blocks
    X <- rpoispp(function(x,y) { 100 * x }, 100, square(1))
    lambda <- 100 * X$x
    Kin <- Kinhom(X, lambda, correction=c("none", "border"))
    lambda2 <- outer(lambda, lambda, "*")
    Ki2 <- Kinhom(X, lambda2=lambda2, diagonal=FALSE,
                  correction=c("translate", "isotropic"))
  }
  if(ALWAYS) {
    #' edge corrections
    rr <- rep(0.1, npoints(Cells))
    eC <- edge.Ripley(Cells, rr)
    eI <- edge.Ripley(Cells, rr, method="interpreted")
    if(max(abs(eC-eI)) > 0.1)
      stop("Ripley edge correction results do not match")
  }
  if(FULLTEST) {
    a <- rmax.Ripley(square(1))
    a <- rmax.Rigid(square(1))
    a <- rmax.Ripley(as.polygonal(square(1)))
    a <- rmax.Rigid(as.polygonal(square(1)))
    a <- rmax.Ripley(letterR)
    a <- rmax.Rigid(letterR)
  }
  if(ALWAYS) {
    #' run slow code for edge correction and compare results
    op <- spatstat.options(npixel=128)
    X <- Redwood[c(TRUE, FALSE, FALSE, FALSE)]
    Window(X) <- as.polygonal(Window(X))
    Eapprox <- edge.Trans(X)
    Eexact <- edge.Trans(X, exact=TRUE)
    maxrelerr <- max(abs(1 - range(Eapprox/Eexact)))
    if(maxrelerr > 0.1)
      stop(paste("Exact and approximate algorithms for edge.Trans disagree by",
                 paste0(round(100*maxrelerr), "%")),
           call.=FALSE)
    spatstat.options(op)
  }
})

local({
  if(FULLTEST) {
    #' ----  multitype ------
    K <- Kcross(Amacrine, correction=c("none", "bord.modif"))
    K <- Kcross(Amacrine, correction=c("none", "bord", "bord.modif"),
                          ratio=TRUE)
    #' inhomogeneous multitype
    K2 <- Kcross.inhom(Amacrine, lambdaX=densityfun(Amacrine))
    K3 <- Kcross.inhom(Amacrine, lambdaX=density(Amacrine, at="points"))
    K5 <- Kcross.inhom(Amacrine, correction="bord.modif")
    #' markconnect, markcorr
    M <- markconnect(Amacrine, "on", "off", normalise=TRUE)
    M <- markcorr(longleaf, normalise=TRUE,
                  correction=c("isotropic", "translate", "border", "none"))
    M <- markcorr(longleaf, normalise=TRUE, fargs=list())
    #' Kmark (=markcorrint)
    X <- runifpoint(100) %mark% runif(100)
    km <- Kmark(X, f=atan2)
    km <- Kmark(X, f1=sin)
    km <- Kmark(X, f="myfun")
    aa <- Kmark(X, normalise=FALSE, returnL=FALSE)
    aa <- Kmark(X, normalise=FALSE, returnL=TRUE)
    aa <- Kmark(X, normalise=TRUE,  returnL=FALSE)
    aa <- Kmark(X, normalise=TRUE,  returnL=TRUE)
  }
})

local({
  if(FULLTEST) {
    #'    various modified K functions
    #'
    #'   directional K functions
    #'
    a <- Ksector(swedishpines,
                 -pi/2, pi/2, units="radians",
                 correction=c("none", "border", "bord.modif",
                              "Ripley", "translate"),
                 ratio=TRUE)
    plot(a)
    #'
    #'   local K functions
    #'
    Z <- as.im(intensity(swedishpines), W=Window(swedishpines))
    ZX <- Z[swedishpines]
    a <- localLinhom(swedishpines, lambda=Z)
    a <- localLinhom(swedishpines, lambda=ZX)
    a <- localLinhom(swedishpines, lambda=Z, correction="none")
    a <- localLinhom(swedishpines, lambda=Z, correction="translate")
    a <- localLcross(Amacrine)
    a <- localLcross(Amacrine, from="off", to="off")
    a <- localKdot(Amacrine)
    a <- localLdot(Amacrine)
    a <- localKcross.inhom(Amacrine)
    a <- localLcross.inhom(Amacrine)
    Zed <- solapply(intensity(amacrine), as.im, W=Window(amacrine))
    Lum <- evaluateCovariateAtPoints(Zed, Amacrine)
    moff <- (marks(Amacrine) == "off")
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed)
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Lum)
    a <- localLcross.inhom(Amacrine, from="off", to="on",
                           lambdaFrom=Lum[moff], lambdaTo=Lum[!moff])
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed,
                           correction="none")
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed,
                           correction="translate")
    #'
    #' cases of resolve.lambdacross
    #'
    h <- resolve.lambdacross(Amacrine, moff, !moff)
    h <- resolve.lambdacross(Amacrine, moff, !moff, lambdaX=Zed)
    h <- resolve.lambdacross(Amacrine, moff, !moff, lambdaX=Lum)
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaI=Zed[["off"]], lambdaJ=Zed[["on"]])
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaI=Lum[moff], lambdaJ=Lum[!moff])
    d <- densityfun(unmark(Amacrine), sigma=0.1)
    dm <- lapply(split(Amacrine), densityfun, sigma=0.1)
    h <- resolve.lambdacross(Amacrine, moff, !moff, lambdaX=d)
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaI=dm[["off"]], lambdaJ=dm[["on"]])
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaX=function(x,y,m){ d(x,y) })
    #'
    #' multitype inhomogeneous pcf
    #'
    g <- pcfcross.inhom(Amacrine, 
                        lambdaI=dm[["off"]], lambdaJ=dm[["on"]])

    #'
    #'   lohboot code blocks
    #'
    Ared <- lohboot(Redwood, fun="Kest", block=TRUE,
                    Vcorrection=TRUE, global=FALSE, correction="none")
    Bred <- lohboot(Redwood, block=TRUE, basicboot=TRUE, global=FALSE)
    Cred <- lohboot(Redwood, fun=Kest, block=TRUE, global=TRUE,
                    correction="translate")
    Dred <- lohboot(Redwood, Lest)
    Kred <- lohboot(Redwood, Kinhom)
    Lred <- lohboot(Redwood, Linhom)
    gred <- lohboot(Redwood, pcfinhom, sigma=0.1)
    #'
    X <- runifpoint(100, letterR)
    AX <- lohboot(X, block=TRUE, nx=7, ny=10)
    #'    multitype
    b <- lohboot(Amacrine, Kcross)
    b <- lohboot(Amacrine, Lcross)
    b <- lohboot(Amacrine, Kdot)
    b <- lohboot(Amacrine, Ldot)
    b <- lohboot(Amacrine, Kcross.inhom)
    b <- lohboot(Amacrine, Lcross.inhom)
    
    ##  Kscaled
    A <- Lscaled(japanesepines, renormalise=TRUE, correction="all")
  }
})
  
local({
  if(ALWAYS) {
    #' From Ege, in response to a stackoverflow question.
    #' The following example has two points separated by r = 1 with 1/4 of the
    #' circumference outside the 10x10 window (i.e. area 100).
    #' Thus the value of K^(r) should jump from 0 to 
    #' 100/(2\cdot 1)\cdot ((3/4)^{-1} + (3/4)^{-1}) = 100 \cdot 4/3 = 133.333.
    x <- c(4.5,5.5)
    y <- c(10,10)-sqrt(2)/2
    W <- square(10)
    X <- ppp(x, y, W)
    compere <- function(a, b, where, tol=1e-6) {
      descrip <- paste("discrepancy in isotropic edge correction", where)
      err <- as.numeric(a) - as.numeric(b)
      maxerr <- max(abs(err))
      blurb <- paste(descrip, "is", paste0(signif(maxerr, 4), ","), 
                     if(maxerr > tol) "exceeding" else "within",
                     "tolerance of", tol)
      message(blurb)
      if(maxerr > tol) {
        message(paste("Discrepancies:", paste(err, collapse=", ")))
        stop(paste("excessive", descrip), call.=FALSE)
      }
      invisible(TRUE)
    }
    ## Testing:
    eX <- edge.Ripley(X, c(1,1))
    compere(eX, c(4/3,4/3), "at interior point of rectangle")
    ## Corner case:
    Y <- X
    Y$x <- X$x-4.5+sqrt(2)/2
    eY <- edge.Ripley(Y, c(1,1))
    compere(eY, c(2,4/3), "near corner of rectangle")
    ## Invoke polygonal code
    Z <- rotate(Y, pi/4)
    eZdebug <- edge.Ripley(Z, c(1,1), internal=list(debug=TRUE))
    compere(eZdebug, c(2,4/3), "at interior point of polygon (debug on)")
    ## test validity without debugger,in case of quirks of compiler optimisation
    eZ <- edge.Ripley(Z, c(1,1))
    compere(eZ,      c(2,4/3), "at interior point of polygon (debug off)")
  }
})



reset.spatstat.options()
