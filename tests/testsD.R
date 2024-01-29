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
#'  tests/density.R
#'
#'  Test behaviour of density() methods,
#'                    relrisk(), Smooth()
#'                    and inhomogeneous summary functions
#'                    and idw, adaptive.density, intensity
#'                    and SpatialMedian, SpatialQuantile
#'
#'  $Revision: 1.67 $  $Date: 2024/01/29 07:07:16 $
#'

if(!FULLTEST)
  spatstat.options(npixel=32, ndummy.min=16)


local({

  # test all cases of density.ppp and densityfun.ppp
  
  tryit <- function(..., do.fun=TRUE, badones=FALSE) {
    Z <- density(cells, ..., at="pixels")
    Z <- density(cells, ..., at="points")
    if(do.fun) {
      f <- densityfun(cells, ...)
      U <- f(0.1, 0.3)
      if(badones) {
        U2 <- f(1.1, 0.3)
        U3 <- f(1.1, 0.3, drop=FALSE)
      }
    }
    return(invisible(NULL))
  }

  if(ALWAYS) {
    tryit(0.05)
    tryit(0.05, diggle=TRUE)
    tryit(0.05, se=TRUE)
    tryit(0.05, weights=expression(x))
    tryit(0.07, kernel="epa")
    tryit(sigma=Inf)
    tryit(0.05, badones=TRUE)
  }
  if(FULLTEST) {
    tryit(0.07, kernel="quartic")
    tryit(0.07, kernel="disc")
    tryit(0.07, kernel="epa", weights=expression(x))
    tryit(sigma=Inf, weights=expression(x))
  }
  
  V <- diag(c(0.05^2, 0.07^2))

  if(ALWAYS) {
    tryit(varcov=V)
  }
  if(FULLTEST) {
    tryit(varcov=V, diggle=TRUE)
    tryit(varcov=V, weights=expression(x))
    tryit(varcov=V, weights=expression(x), diggle=TRUE)
    Z <- distmap(runifpoint(5, Window(cells)))
    tryit(0.05, weights=Z)
    tryit(0.05, weights=Z, diggle=TRUE)
  }

  trymost <- function(...) tryit(..., do.fun=FALSE) 
  wdf <- data.frame(a=1:42,b=42:1)
  if(ALWAYS) {
    trymost(0.05, weights=wdf)
    trymost(sigma=Inf, weights=wdf)
  }
  if(FULLTEST) {
    trymost(0.05, weights=wdf, diggle=TRUE)
    trymost(varcov=V, weights=wdf)
    trymost(varcov=V, weights=expression(cbind(x,y)))
  }

  ## check conservation of mass
  checkconserve <- function(X, xname, sigma, toler=0.01) {
    veritas <- npoints(X)
    vino <- integral(density(X, sigma, diggle=TRUE))
    relerr <- abs(vino - veritas)/veritas
    if(relerr > toler)
      stop(paste("density.ppp(diggle=TRUE) fails to conserve mass:",
                 vino, "!=", veritas,
                 "for", sQuote(xname)),
           call.=FALSE)
    return(relerr)
  }
  if(FULLTEST) {
    checkconserve(cells, "cells", 0.15)
  }
  if(ALWAYS) {
    checkconserve(split(chorley)[["lung"]], "lung", 2)
  }
  
  ## run C algorithm 'denspt'
  opa <- spatstat.options(densityC=TRUE, densityTransform=FALSE)
  if(ALWAYS) {
    tryit(varcov=V)
  }
  if(FULLTEST) {
    tryit(varcov=V, weights=expression(x))
    trymost(varcov=V, weights=wdf)
  }
  spatstat.options(opa)

  crossit <- function(..., sigma=NULL) {
    U <- runifpoint(20, Window(cells))
    a <- densitycrossEngine(cells, U, ..., sigma=sigma)
    a <- densitycrossEngine(cells, U, ..., sigma=sigma, diggle=TRUE)
    invisible(NULL)
  }
  if(ALWAYS) {
    crossit(varcov=V, weights=cells$x)
    crossit(sigma=Inf)
  }
  if(FULLTEST) {
    crossit(varcov=V, weights=wdf)
    crossit(sigma=0.1, weights=wdf)
    crossit(sigma=0.1, kernel="epa", weights=wdf)
  }
  
  ## apply different discretisation rules
  if(ALWAYS) {
    Z <- density(cells, 0.05, fractional=TRUE)
  }
  if(FULLTEST) {
    Z <- density(cells, 0.05, preserve=TRUE)
    Z <- density(cells, 0.05, fractional=TRUE, preserve=TRUE)
  }
        
  ## compare results with different algorithms
  crosscheque <- function(expr) {
    e <- as.expression(substitute(expr))
    ename <- sQuote(deparse(substitute(expr)))
    ## interpreted R
    opa <- spatstat.options(densityC=FALSE, densityTransform=FALSE)
    val.interpreted <- eval(e)
    ## established C algorithm 'denspt'
    spatstat.options(densityC=TRUE, densityTransform=FALSE)
    val.C <- eval(e)
    ## new C algorithm 'Gdenspt' using transformed coordinates
    spatstat.options(densityC=TRUE, densityTransform=TRUE)
    val.Transform <- eval(e)
    spatstat.options(opa)
    if(max(abs(val.interpreted - val.C)) > 0.001)
      stop(paste("Numerical discrepancy between R and C algorithms in",
                 ename))
    if(max(abs(val.C - val.Transform)) > 0.001)
      stop(paste("Numerical discrepancy between C algorithms",
                 "using transformed and untransformed coordinates in",
                 ename))
    invisible(NULL)
  }

  ## execute & compare results of density(at="points") with different algorithms
  wdfr <- cbind(1:npoints(redwood), 2)
  if(ALWAYS) {
    crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE))
    crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE,
                        weights=wdfr[,1]))
    crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE,
                        weights=wdfr))
  }

  ## correctness of non-Gaussian kernel calculation
  leavein <- function(ker, maxd=0.025) {
    ZI <- density(redwood, 0.12, kernel=ker, edge=FALSE,
                  dimyx=256)[redwood]
    ZP <- density(redwood, 0.12, kernel=ker, edge=FALSE,
                  at="points", leaveoneout=FALSE)
    discrep <- max(abs(ZP - ZI))/npoints(redwood)
    if(discrep > maxd) 
      stop(paste("Discrepancy",
                 signif(discrep, 3),
                 "in calculation for", ker, "kernel"))
    return(invisible(NULL))
  }
  if(ALWAYS) {
    leavein("epanechnikov", 0.015)
  }
  if(FULLTEST) {
    leavein("quartic",      0.010)
    leavein("disc",         0.100)
  }

  ## bandwidth selection code blocks
  sigvec <- 0.01 * 2:15
  sigran <- range(sigvec)
  if(ALWAYS) {
    bw.ppl(redwood, sigma=sigvec)
    bw.CvL(redwood, sigma=sigvec)
  }
  if(FULLTEST) {
    bw.ppl(redwood, srange=sigran, ns=5)
    bw.CvL(redwood, srange=sigran, ns=5)
  }
  ## adaptive bandwidth
  if(ALWAYS) {
    a <- bw.abram(redwood)
  }
  if(FULLTEST) {
    a <- bw.abram(redwood, pilot=density(redwood, 0.2))
    a <- bw.abram(redwood, smoother="densityVoronoi", at="pixels")
  }
  
  ## Kinhom
  if(ALWAYS) {
    lam <- density(redwood)
    K <- Kinhom(redwood, lam)
  
    lamX <- density(redwood, at="points")
    KX <- Kinhom(redwood, lamX)
  }

  ## test all code cases of new 'relrisk.ppp' algorithm
  pants <- function(..., X=ants, sigma=100, se=TRUE) {
    a <- relrisk(X, sigma=sigma, se=se, ...)
    return(TRUE)
  }
  if(ALWAYS) {
    pants()
    pants(diggle=TRUE)
    pants(edge=FALSE)
    pants(at="points")
    pants(casecontrol=FALSE)
    pants(relative=TRUE)
    pants(sigma=Inf)
    pants(sigma=NULL, varcov=diag(c(100,100)^2))
    f <- 1/area(Window(ants))
    pants(fudge=f)
  }
  if(FULLTEST) {
    pants(diggle=TRUE, at="points")
    pants(edge=FALSE, at="points", fudge=f)
    pants(casecontrol=FALSE, relative=TRUE)
    pants(casecontrol=FALSE,at="points")
    pants(relative=TRUE,at="points", fudge=f)
    pants(casecontrol=FALSE, relative=TRUE,at="points")
    pants(relative=TRUE, control="Cataglyphis", case="Messor", fudge=f)
    pants(relative=TRUE, control="Cataglyphis", case="Messor", at="points")
    pants(casecontrol=FALSE, case="Messor", se=FALSE)
    pants(case=2, at="pixels", relative=TRUE)
    pants(case=2, at="points", relative=TRUE)
    pants(case=2, at="pixels", relative=FALSE)
    pants(case=2, at="points", relative=FALSE)
  }
  if(ALWAYS) {
    ## underflow example from stackoverflow!
    funky <- scanpp("funky.tab", owin(c(4, 38), c(0.3, 17)))
    P <- relrisk(funky, 0.5)
    R <- relrisk(funky, 0.5, relative=TRUE)
  }
  ## more than 2 types
  if(ALWAYS) {
    pants(X=sporophores)
    pants(X=sporophores, sigma=20, at="points")
    pants(X=sporophores, sigma=20, at="points", fudge=f)
    bw.relrisk(sporophores, method="leastsquares")
  }
  if(FULLTEST) {
    pants(X=sporophores, sigma=20, relative=TRUE, at="points", fudge=f)
    pants(X=sporophores, sigma=20, at="pixels", se=FALSE)
    pants(X=sporophores, sigma=20, relative=TRUE, at="pixels", se=FALSE)
    bw.relrisk(sporophores, method="weightedleastsquares")
  }
  
  ## execute Smooth.ppp and Smoothfun.ppp in all cases
  stroke <- function(..., Y = longleaf, FUN=TRUE) {
    Z <- Smooth(Y, ..., at="pixels")
    Z <- Smooth(Y, ..., at="points", leaveoneout=TRUE)
    Z <- Smooth(Y, ..., at="points", leaveoneout=FALSE)
    if(FUN) {
      f <- Smoothfun(Y, ...)
      f(120, 80)
      f(Y[1:2])
      f(Y[FALSE])
      U <- as.im(f)
    }
    return(invisible(NULL))
  }
  if(ALWAYS) {
    stroke()
    stroke(5, diggle=TRUE)
    stroke(5, geometric=TRUE)
    stroke(1e-6) # generates warning about small bandwidth
    stroke(5, weights=expression(x))
    stroke(5, kernel="epa")
    stroke(sigma=Inf)
    stroke(varcov1=diag(c(1,1))) # 'anisotropic' code
  }
  if(FULLTEST) {
    Z <- as.im(function(x,y){abs(x)+1}, Window(longleaf))
    stroke(5, weights=Z)
    stroke(5, weights=runif(npoints(longleaf)))
    stroke(varcov=diag(c(25, 36)))
    stroke(varcov=diag(c(25, 36)), weights=runif(npoints(longleaf)))
    stroke(5, Y=longleaf %mark% 1)
    stroke(5, Y=cut(longleaf,breaks=3))
    stroke(5, weights=Z, geometric=TRUE)
    g <- function(x,y) { dnorm(x, sd=10) * dnorm(y, sd=10) }
    stroke(kernel=g, cutoff=30, FUN=FALSE)
    stroke(kernel=g, cutoff=30, scalekernel=TRUE, sigma=1, FUN=FALSE)
  }
  if(FULLTEST) {
    ## standard errors - single column of marks
    stroke(sigma=5, se=TRUE)
    stroke(sigma=5, se=TRUE, loctype="f")
    w <- runif(npoints(longleaf))
    stroke(sigma=5, se=TRUE, weights=w, loctype="r", wtype="i")
    stroke(sigma=5, se=TRUE, weights=w, loctype="r", wtype="m")
    stroke(sigma=5, se=TRUE, weights=w, loctype="f", wtype="i")
    stroke(sigma=5, se=TRUE, weights=w, loctype="f", wtype="m")
  }
  
  niets <- markmean(longleaf, 9)
  
  strike <- function(..., Y=finpines) {
    Z <- Smooth(Y, ..., at="pixels")
    Z <- Smooth(Y, ..., at="points", leaveoneout=TRUE)
    Z <- Smooth(Y, ..., at="points", leaveoneout=FALSE)
    f <- Smoothfun(Y, ...)
    f(4, 1)
    f(Y[1:2])
    f(Y[FALSE])
    U <- as.im(f)
    return(invisible(NULL))
  }
  if(ALWAYS) {
    strike()
    strike(sigma=1.5, kernel="epa")
    strike(varcov=diag(c(1.2, 2.1)))
    strike(sigma=1e-6)
    strike(sigma=Inf)
  }
  if(FULLTEST) {
    strike(sigma=1e-6, kernel="epa")
    strike(1.5, weights=runif(npoints(finpines)))
    strike(1.5, weights=expression(y))
    strike(1.5, geometric=TRUE)
    strike(1.5, Y=finpines[FALSE])
    flatfin <- finpines %mark% data.frame(a=rep(1, npoints(finpines)), b=2)
    strike(1.5, Y=flatfin)
    strike(1.5, Y=flatfin, geometric=TRUE)
  }
  if(FULLTEST) {
    ## standard errors - multivariate marks
    strike(sigma=1.5, se=TRUE)
    strike(sigma=1.5, se=TRUE, loctype="f")
    w <- runif(npoints(finpines))
    strike(sigma=1.5, se=TRUE, weights=w, loctype="r", wtype="i")
    strike(sigma=1.5, se=TRUE, weights=w, loctype="r", wtype="m")
    strike(sigma=1.5, se=TRUE, weights=w, loctype="f", wtype="i")
    strike(sigma=1.5, se=TRUE, weights=w, loctype="f", wtype="m")
  }
  opx <- spatstat.options(densityTransform=FALSE)
  if(ALWAYS) {
    stroke(5, Y=longleaf[order(longleaf$x)], sorted=TRUE)
  }
  if(FULLTEST) {
    strike(1.5, Y=finpines[order(finpines$x)], sorted=TRUE)
  }
  spatstat.options(opx)

  ## detect special cases
  if(ALWAYS) {
    Smooth(longleaf[FALSE])
    Smooth(longleaf, minnndist(longleaf))
    Xconst <- cells %mark% 1
    Smooth(Xconst, 0.1)
    Smooth(Xconst, 0.1, at="points")
    Smooth(cells %mark% runif(42), sigma=Inf)
    Smooth(cells %mark% runif(42), sigma=Inf, at="points")
    Smooth(cells %mark% runif(42), sigma=Inf, at="points", leaveoneout=FALSE)
    Smooth(cut(longleaf, breaks=4))
  }
  
  ## code not otherwise reached
  if(ALWAYS) {
    smoothpointsEngine(cells, values=rep(1, npoints(cells)), sigma=0.2)
  }
  if(FULLTEST) {
    smoothpointsEngine(cells, values=runif(npoints(cells)), sigma=Inf)
    smoothpointsEngine(cells, values=runif(npoints(cells)), sigma=1e-16)
  }
  
  ## validity of Smooth.ppp(at='points')
  Y <- longleaf %mark% runif(npoints(longleaf), min=41, max=43)
  Z <- Smooth(Y, 5, at="points", leaveoneout=TRUE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=TRUE)")

  Z <- Smooth(Y, 5, at="points", leaveoneout=FALSE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=FALSE)")

  ## compare Smooth.ppp results with different algorithms
  if(ALWAYS) {
    crosscheque(Smooth(longleaf, at="points", sigma=6))
    wt <- runif(npoints(longleaf))
    crosscheque(Smooth(longleaf, at="points", sigma=6, weights=wt))
  }
  if(FULLTEST) {
    vc <- diag(c(25,36))
    crosscheque(Smooth(longleaf, at="points", varcov=vc))
    crosscheque(Smooth(longleaf, at="points", varcov=vc, weights=wt))
  }
  ## drop-dimension coding errors
  if(FULLTEST) {
    X <- longleaf
    marks(X) <- cbind(marks(X), 1)
    Z <- Smooth(X, 5)

    ZZ <- bw.smoothppp(finpines, hmin=0.01, hmax=0.012, nh=2) # reshaping problem
  }

  ## geometric-mean smoothing
  if(ALWAYS) {
    U <- Smooth(longleaf, 5, geometric=TRUE)
  }
  if(FULLTEST) {
    UU <- Smooth(X, 5, geometric=TRUE)
    V <- Smooth(longleaf, 5, geometric=TRUE, at="points")
    VV <- Smooth(X, 5, geometric=TRUE, at="points")
  }

  if(FULLTEST) {
    ## isotropic and anisotropic cases of bw.smoothppp
    bi <- bw.smoothppp(longleaf)
    ba <- bw.smoothppp(longleaf, varcov1=diag(c(1,1)))
    ## should be equal
    if(abs(bi-ba) > 0.001)
      stop(paste("Inconsistency in bw.smoothppp: isotropic =", bi,
                 "!=", ba, "= anisotropic"))
  }
})

reset.spatstat.options()

local({
  if(ALWAYS) {
    #' Kmeasure, second.moment.engine
    #' Expansion of window
    Zno  <- Kmeasure(redwood, sigma=0.2, expand=FALSE)
    Zyes <- Kmeasure(redwood, sigma=0.2, expand=TRUE)
    #' All code blocks
    sigmadouble <- rep(0.1, 2)
    diagmat <- diag(sigmadouble^2)
    generalmat <- matrix(c(1, 0.5, 0.5, 1)/100, 2, 2)
    Z <- Kmeasure(redwood, sigma=sigmadouble)
    Z <- Kmeasure(redwood, varcov=diagmat)
    Z <- Kmeasure(redwood, varcov=generalmat)
    A <- second.moment.calc(redwood, 0.1, what="all", debug=TRUE)
    B <- second.moment.calc(redwood, varcov=diagmat,    what="all")
    B <- second.moment.calc(redwood, varcov=diagmat,    what="all")
    D <- second.moment.calc(redwood, varcov=generalmat, what="all")
    PR <- pixellate(redwood)
    DRno  <- second.moment.calc(PR, 0.2, debug=TRUE, expand=FALSE,
                                npts=npoints(redwood), obswin=Window(redwood))
    DRyes <- second.moment.calc(PR, 0.2, debug=TRUE, expand=TRUE,
                                npts=npoints(redwood), obswin=Window(redwood))
    DR2 <- second.moment.calc(solist(PR, PR), 0.2, debug=TRUE, expand=TRUE,
                              npts=npoints(redwood), obswin=Window(redwood))
    Gmat <- generalmat * 100
    isoGauss <- function(x,y) {dnorm(x) * dnorm(y)}
    ee <- evaluate2Dkernel(isoGauss, runif(10), runif(10),
                           varcov=Gmat, scalekernel=TRUE)
    isoGaussIm <- as.im(isoGauss, square(c(-3,3)))
    gg <- evaluate2Dkernel(isoGaussIm, runif(10), runif(10),
                           varcov=Gmat, scalekernel=TRUE)
    ## experimental code
    op <- spatstat.options(developer=TRUE)
    DR <- density(redwood, 0.1)
    spatstat.options(op)
  }
})

local({
  if(FULLTEST) {
    #' bandwidth selection
    op <- spatstat.options(n.bandwidth=8)
    bw.diggle(cells) 
    bw.diggle(cells, method="interpreted") # undocumented test
    ##  bw.relrisk(urkiola, hmax=20) is tested in man/bw.relrisk.Rd
    bw.relrisk(urkiola, hmax=20, method="leastsquares")
    bw.relrisk(urkiola, hmax=20, method="weightedleastsquares")
    ZX <- density(swedishpines, at="points")
    bw.pcf(swedishpines, lambda=ZX)
    bw.pcf(swedishpines, lambda=ZX,
           bias.correct=FALSE, simple=FALSE, cv.method="leastSQ")
    spatstat.options(op)
  }
})

local({
  if(FULLTEST) {
    ## idw
    Z <- idw(longleaf, power=4)
    Z <- idw(longleaf, power=4, se=TRUE)
    ZX <- idw(longleaf, power=4, at="points")
    ZX <- idw(longleaf, power=4, at="points", se=TRUE)
  }
  if(ALWAYS) {
    ## former bug in densityVoronoi.ppp 
    X <- redwood[1:2]
    A <- densityVoronoi(X, f=0.51, counting=FALSE, fixed=FALSE, nrep=50, verbose=FALSE)
    ## dodgy code blocks in densityVoronoi.R
    A <- adaptive.density(nztrees, nrep=2, f=0.5, counting=TRUE)
    B <- adaptive.density(nztrees, nrep=2, f=0.5, counting=TRUE, fixed=TRUE)
    D <- adaptive.density(nztrees, nrep=2, f=0.5, counting=FALSE)
    E <- adaptive.density(nztrees, nrep=2, f=0.5, counting=FALSE, fixed=TRUE)
  }
  if(FULLTEST) {
    #' adaptive kernel estimation
    d10 <- nndist(nztrees, k=10)
    d10fun <- distfun(nztrees, k=10)
    d10im  <- as.im(d10fun)
    uN <- 2 * runif(npoints(nztrees))
    AA <- densityAdaptiveKernel(nztrees, bw=d10)
    BB <- densityAdaptiveKernel(nztrees, bw=d10, weights=uN)
    DD <- densityAdaptiveKernel(nztrees, bw=d10fun, weights=uN)
    EE <- densityAdaptiveKernel(nztrees, bw=d10im, weights=uN)
  }
})

local({
  if(FULLTEST) {
    ## cases of 'intensity' etc
    a <- intensity(amacrine, weights=expression(x))
    SA <- split(amacrine)
    a <- intensity(SA, weights=expression(x))
    a <- intensity(SA, weights=amacrine$x)

    ## check infrastructure for 'densityfun'
    f <- densityfun(cells, 0.05)
    Z <- as.im(f)
    Z <- as.im(f, W=square(0.5))
  }
})

local({
  if(FULLTEST) {
    ## other cases of SpatialQuantile.ppp
    X <- longleaf
    marks(X) <- round(marks(X), -1)
    Z <- SpatialMedian(X, 30, type=4)
    ZX <- SpatialMedian(X, 30, type=4, at="points")
    ZXP <- SpatialMedian(X, 30, at="points", leaveoneout=FALSE)
  }
})



reset.spatstat.options()

