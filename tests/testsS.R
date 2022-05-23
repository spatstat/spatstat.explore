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
#'   tests/sdr.R
#'
#'   $Revision: 1.2 $ $Date: 2020/05/01 09:59:59 $

if(FULLTEST) {
local({
  AN <- sdr(bei, bei.extra, method="NNIR")
  AV <- sdr(bei, bei.extra, method="SAVE")
  AI <- sdr(bei, bei.extra, method="SIR")
  AT <- sdr(bei, bei.extra, method="TSE")
  subspaceDistance(AN$B, AV$B)
  dimhat(AN$M)
})
}
##
##  tests/segments.R
##   Tests of psp class and related code
##                      [SEE ALSO: tests/xysegment.R]
##
##  $Revision: 1.33 $  $Date: 2022/05/22 08:39:47 $


local({
  if(ALWAYS) { # C code
    #' tests of density.psp
    Y <- edges(letterR)
    Window(Y) <- grow.rectangle(Frame(Y), 0.4)
    YC <- density(Y, 0.2, method="C", edge=FALSE, dimyx=64)
    YI <- density(Y, 0.2, method="interpreted", edge=FALSE, dimyx=64)
    YF <- density(Y, 0.2, method="FFT", edge=FALSE, dimyx=64)
    xCI <- max(abs(YC/YI - 1))
    xFI <- max(abs(YF/YI - 1))
    cat(paste("xCI =", xCI, "\txFI =", signif(xFI, 5)), fill=TRUE)
    if(xCI > 0.01) stop(paste("density.psp C algorithm relative error =", xCI))
    if(xFI > 0.1) stop(paste("density.psp FFT algorithm relative error =", xFI))

    B <- square(0.3)
    density(Y, 0.2, at=B)
    density(Y, 0.2, at=B, edge=TRUE, method="C")
    Z <- runifrect(3, B)
    density(Y, 0.2, at=Z)
    density(Y, 0.2, at=Z, edge=TRUE, method="C")
  }

  if(FULLTEST) {
    #' segment clipping in window (bug found by Rolf)
    set.seed(42)
    X <- runifpoint(50, letterR)
    SP <- dirichletEdges(X) #' clip to polygonal window
    Window(X) <- as.mask(Window(X))
    SM <- dirichletEdges(X) #' clip to mask window
  }
  
  if(FULLTEST) {
    #' test rshift.psp and append.psp with marks (Ute Hahn)
    m <- data.frame(A=1:10, B=letters[1:10])
    g <- gl(3, 3, length=10)
    X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)
    Y <- rshift(X, radius = 0.1)
    Y <- rshift(X, radius = 0.1, group=g)
    #' mark management
    b <- data.frame(A=1:10)
    X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=b)
    stopifnot(is.data.frame(marks(X)))
    Y <- rshift(X, radius = 0.1)
    Y <- rshift(X, radius = 0.1, group=g)
  }

})



#
## tests/sigtraceprogress.R
#
## Tests of *.sigtrace and *.progress
#
## $Revision: 1.5 $ $Date: 2020/05/01 09:59:59 $

if(FULLTEST) {
local({
  plot(dclf.sigtrace(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dclf.progress(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.sigtrace(redwood, nsim=5, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.progress(redwood, nsim=5, alternative="greater", rmin=0.02,
                   verbose=FALSE))
  ## test 'leave-two-out' algorithm
  a <- dclf.sigtrace(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                     verbose=FALSE)
  aa <- dclf.progress(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                      verbose=FALSE)
  b <- dg.sigtrace(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2)
  bb <- dg.progress(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2,
                    verbose=FALSE)
  ## other code blocks
  e <- mad.progress(redwood, nsim=5)
  e <- mad.progress(redwood, nsim=19, alpha=0.05)
  f <- dclf.progress(redwood, nsim=5, scale=function(x) x^2)
  f <- dclf.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
  g <- dg.progress(redwood, nsim=5, scale=function(x) x^2)
  g <- dg.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
})
}
#'
#'    tests/ssf.R
#'
#'   Tests of 'ssf' class
#'
#'   $Revision: 1.5 $ $Date: 2020/12/04 08:02:25 $
#'

if(FULLTEST) {
local({
  Y <- cells[1:5]
  X <- rsyst(Window(Y), 5)
  Z <- runifpoint(3, Window(Y))
  f1 <- ssf(X, nncross(X,Y,what="dist"))
  f2 <- ssf(X, nncross(X,Y,what="dist", k=1:2))
  image(f1)
  g1 <- as.function(f1)
  g1(Z)
  g2 <- as.function(f2)
  g2(Z)
  plot(f1, style="contour")
  plot(f1, style="imagecontour")
  contour(f1)
  apply.ssf(f2, 1, sum)
  range(f1)
  min(f1)
  max(f1)
  integral(f1, weights=tile.areas(dirichlet(X)))
})
}
#'
#'   tests/sumfun.R
#'
#'   Tests of code for summary functions
#'
#'   $Revision: 1.9 $ $Date: 2022/05/22 08:45:23 $

if(ALWAYS) { # involves C code 
local({
  W <- owin(c(0,1), c(-1/2, 0))
  Gr <- Gest(redwood, correction="all",domain=W)
  Fr <- Fest(redwood, correction="all",domain=W)
  Jr <- Jest(redwood, correction="all",domain=W)
  
  F0 <- Fest(redwood[FALSE], correction="all")
  Fh <- Fest(humberside, domain=erosion(Window(humberside), 100))

  FIr <- Finhom(redwood, savelambda=TRUE, ratio=TRUE)
  JIr <- Jinhom(redwood, savelambda=TRUE, ratio=TRUE)
  
  Ga <- Gcross(amacrine, correction="all")
  Ia <- Iest(amacrine, correction="all")
  lam <- intensity(amacrine)
  lmin <- 0.9 * min(lam)
  nJ <- sum(marks(amacrine) == "off")
  FM <- FmultiInhom(amacrine, marks(amacrine) == "off",
                    lambdaJ=rep(lam["off"], nJ),
                    lambdamin = lmin)
  GM <- GmultiInhom(amacrine, marks(amacrine) == "on",
                    marks(amacrine) == "off",
                    lambda=lam[marks(amacrine)],
                    lambdamin=lmin,
                    ReferenceMeasureMarkSetI=42)

  a <- compileCDF(D=nndist(redwood),
                  B=bdist.points(redwood),
                  r=seq(0, 1, length=256))

  #' Tstat (triplet) function, all code blocks
  a <- Tstat(redwood, ratio=TRUE,
             correction=c("none", "border", "bord.modif", "translate"))
  
  ## distance argument spacing and breakpoints
  e <- check.finespacing(c(0,1,2), eps=0.1, action="silent")
  b <- as.breakpts(pi, 20)
  b <- as.breakpts(42, max=pi, npos=20)
  b <- even.breaks.owin(letterR)
})
}
