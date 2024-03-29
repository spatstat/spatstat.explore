#
# bermantest.R
#
# Test statistics from Berman (1986)
#
#  $Revision: 1.26 $  $Date: 2023/06/20 02:43:35 $
#
#


berman.test <- function(...) {
  UseMethod("berman.test")
}

berman.test.ppp <-
  function(X, covariate,
           which=c("Z1", "Z2"),
           alternative=c("two.sided", "less", "greater"),
           ...) {
    Xname <- short.deparse(substitute(X))
    covname <- short.deparse(substitute(covariate))
    force(covariate)
    if(is.character(covariate)) covname <- covariate
    which <- match.arg(which)
    alternative <- match.arg(alternative)

    fitcsr <- exactppm(X)
    dont.complain.about(fitcsr)
    
    do.call(bermantestEngine,
            resolve.defaults(list(quote(fitcsr),
                                  quote(covariate),
                                  which, alternative),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

## Code for berman.test.ppm is moved to spatstat.model


bermantestEngine <- function(model, covariate,
                             which=c("Z1", "Z2"),
                             alternative=c("two.sided", "less", "greater"),
                             ...,
                             modelname, covname, dataname="") {

  csr <- is.poisson(model) && is.stationary(model)
  if(missing(modelname))
    modelname <- if(csr) "CSR" else short.deparse(substitute(model))
  if(missing(covname)) {
    covname <- short.deparse(substitute(covariate))
    if(is.character(covariate)) covname <- covariate
  }

  which <- match.arg(which)
  alternative <- match.arg(alternative)

  if(!is.poisson(model))
    stop("Only implemented for Poisson point process models")

  #'  compute required data 
  fram <- spatialCDFframe(model, covariate, ...,
                          modelname=modelname,
                          covname=covname,
                          dataname=dataname)
  #'  evaluate berman test statistic 
  result <- bermantestCalc(fram, which=which, alternative=alternative)

  return(result)
}

bermantestCalc <- function(fram,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {

  which <- match.arg(which)
  alternative <- match.arg(alternative)

  verifyclass(fram, "spatialCDFframe")
  fvalues <- fram$values
  info    <- fram$info
  
  ## values of covariate at data points
  ZX <- fvalues$ZX
  ## transformed to Unif[0,1] under H0
  U  <- fvalues$U
  ## values of covariate at pixels
  Zvalues <- fvalues$Zvalues
  ## corresponding pixel areas/weights
  weights <- fvalues$weights
  ## intensity of model
  lambda  <- fvalues$lambda

  ## names 
  modelname <- info$modelname
  dataname  <- info$dataname
  covname   <- info$covname
  
  switch(which,
         Z1={
           #......... Berman Z1 statistic .....................
           method <-
             paste("Berman Z1 test of",
                   if(info$csr) "CSR" else "inhomogeneous Poisson process",
                   "in", info$spacename)
           # sum of covariate values at data points
           Sn <- sum(ZX)
           # predicted mean and variance
           lamwt <- lambda * weights
           En    <- sum(lamwt)
           ESn   <- sum(lamwt * Zvalues)
           varSn <- sum(lamwt * Zvalues^2)
           # working, for plot method
           working <- list(meanZX=mean(ZX),
                           meanZ=ESn/En)
           # standardise
           statistic <- (Sn - ESn)/sqrt(varSn)
           names(statistic) <- "Z1"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="mean value of covariate at random points is less than predicted under model",
                              greater="mean value of covariate at random points is greater than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname))
         },
         Z2={
           #......... Berman Z2 statistic .....................
           method <-
             paste("Berman Z2 test of",
                   if(info$csr) "CSR" else "inhomogeneous Poisson process",
                   "in", info$spacename)
           npts <- length(ZX)
           statistic <- sqrt(12/npts) * (sum(U) - npts/2)
           working <- list(meanU=mean(U))
           names(statistic) <- "Z2"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="covariate values at random points have lower quantiles than predicted under model",
                              greater="covariate values at random points have higher quantiles than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname), "\n\t",
                              "and transformed to uniform distribution under",
                              if(info$csr) modelname else sQuote(modelname))
         })
           
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=method,
              which=which,
              working=working,
              data.name=valuename,
              fram=fram)
  class(out) <- c("htest", "bermantest")
  return(out)
}

plot.bermantest <-
  function(x, ...,
           lwd=par("lwd"), col=par("col"), lty=par("lty"),
           lwd0=lwd, col0=2, lty0=2)
{
  fram <- x$fram
  if(!is.null(fram)) {
    values <- fram$values
    info <- fram$info
  } else {
    # old style
    ks <- x$ks
    values <- attr(ks, "prep")
    info <- attr(ks, "info")
  }
  work <- x$working
  op <- options(useFancyQuotes=FALSE)
  on.exit(options(op))
  switch(x$which,
         Z1={
           # plot cdf's of Z
           FZ <- values$FZ
           xxx <- get("x", environment(FZ))
           yyy <- get("y", environment(FZ))
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z1 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call(plot.default,
                   resolve.defaults(
                                    list(x=xxx, y=yyy, type="l"),
                                    list(...),
                                    list(lwd=lwd0, col=col0, lty=lty0),
                                    list(xlab=info$covname,
                                         ylab="probability",
                                         main=main)))
           FZX <- values$FZX
           if(is.null(FZX))
             FZX <- ecdf(values$ZX)
           plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
           abline(v=work$meanZ, lwd=lwd0,col=col0, lty=lty0, xpd=FALSE)
           abline(v=work$meanZX, lwd=lwd,col=col, lty=lty, xpd=FALSE)
         },
         Z2={
           # plot cdf of U
           U <- values$U
           cdfU <- ecdf(U)
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z2 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           dont.complain.about(cdfU)
           do.call(plot.ecdf,
                   resolve.defaults(
                                    list(quote(cdfU)),
                                    list(...),
                                    list(do.points=FALSE, asp=1),
                                    list(xlim=c(0,1), ylim=c(0,1), pty="s"),
                                    list(lwd=lwd, col=col, lty=lty),
                                    list(xlab="U", ylab="relative frequency"),
                                    list(main=main)))
           abline(0,1,lwd=lwd0,col=col0,lty=lty0, xpd=FALSE)
           abline(v=0.5, lwd=lwd0,col=col0,lty=lty0, xpd=FALSE)
           abline(v=work$meanU, lwd=lwd,col=col,lty=lty, xpd=FALSE)
         })
  options(op)
  return(invisible(NULL))
}

