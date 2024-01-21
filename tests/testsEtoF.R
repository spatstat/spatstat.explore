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
#  tests/envelopes.R
#
#  Test validity of envelope data
#
#  $Revision: 1.29 $  $Date: 2024/01/10 13:45:29 $
#

local({
  


## check envelope calls from 'alltypes'
if(ALWAYS) a <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE)
if(FULLTEST) b <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE, global=TRUE)
## check 'transform' idioms
if(ALWAYS) A <- envelope(cells, Kest, nsim=4, transform=expression(. - .x))
if(FULLTEST) B <- envelope(cells, Kest, nsim=4, transform=expression(sqrt(./pi) - .x))


# check conditional simulation
if(FULLTEST) {
  e1 <- envelope(cells, Kest, nsim=4, fix.n=TRUE)
  e2 <- envelope(amacrine, Kest, nsim=4, fix.n=TRUE)
  e3 <- envelope(amacrine, Kcross, nsim=4, fix.marks=TRUE)
  e4 <- envelope(finpines, Kest, nsim=4, fix.n=TRUE) # multiple columns of marks
  e5 <- envelope(finpines, Kest, nsim=4, fix.marks=TRUE)
}


## check pooling of envelopes in global case
E1 <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE)
E2 <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE)
p12 <- pool(E1, E2)
p12 <- pool(E1, E2, savefuns=TRUE)
if(FULLTEST) {
  F1 <- envelope(cells, Kest, nsim=5,
                 savefuns=TRUE, savepatterns=TRUE, global=TRUE)
  F2 <- envelope(cells, Kest, nsim=12,
                 savefuns=TRUE, savepatterns=TRUE, global=TRUE)
  p12 <- pool(F1, F2)
  p12 <- pool(F1, F2, savefuns=TRUE, savepatterns=TRUE)
  E1r <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE,
                  ginterval=c(0.05, 0.15))
  E2r <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE,
                  ginterval=c(0.05, 0.15))
  p12r <- pool(E1r, E2r)
}

if(FULLTEST) {
  #' as.data.frame.envelope
  Nsim <- 5
  E <- envelope(cells, nsim=Nsim, savefuns=TRUE)
  A <- as.data.frame(E)
  B <- as.data.frame(E, simfuns=TRUE)
  stopifnot(ncol(B) - ncol(A) == Nsim)
}

if(FULLTEST) {
  #' cases not covered elsewhere
  A <- envelope(cells, nsim=5, alternative="less",
                do.pwrong=TRUE, use.theory=FALSE,
                savepatterns=TRUE, savefuns=TRUE)
  print(A)
  B <- envelope(A, nsim=5, savefuns=TRUE)
  D <- envelope(cells, "Lest", nsim=5)
  
  UU <- envelope(cells, nsim=5, foreignclass="ppp", clipdata=TRUE)
  
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", global=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="less", global=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", VARIANCE=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", VARIANCE=TRUE)

  #' spotted by Art Stock - bugs in ratfv class support
  BB <- envelope(redwood, Kinhom, nsim=5, sigma=bw.scott,              ratio=TRUE, correction="border")
  CC <- envelope(redwood, Kinhom, nsim=5, sigma=bw.scott, global=TRUE, ratio=TRUE, correction="border")
  DD <- envelope(redwood, Finhom, nsim=5, sigma=bw.scott,              ratio=TRUE, correction="trans")
  EE <- envelope(redwood, Finhom, nsim=5, sigma=bw.scott, global=TRUE, ratio=TRUE, correction="trans")
  
  #'  envelopes based on sample variance
  E <- envelope(cells, nsim=8, VARIANCE=TRUE)
  G <- envelope(cells, nsim=8, VARIANCE=TRUE,
                use.theory=FALSE, do.pwrong=TRUE)
  print(G)
  #' summary method
  summary(E)
  summary(envelope(cells, nsim=5, simulate=expression(runifpoint(42))))
  #' weights argument
  H1 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  H2 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  J1 <- envelope(cells, nsim=4, weights=npoints, VARIANCE=TRUE)
  J2 <- envelope(cells, nsim=4, weights=npoints, VARIANCE=TRUE)
  #' pooling with weights
  H <- pool(H1, H2)
  J <- pool(J1, J2)
  #' pooling envelopes with non-identical attributes
  H0 <- envelope(cells, nsim=4, savefuns=TRUE)
  HH <- pool(H0, H1)
  #' malformed argument 'simulate'
  A <- replicate(3, list(list(runifpoint(ex=cells))))   # list(list(ppp), list(ppp), list(ppp))
  E <- envelope(cells, simulate=A, nsim=3)
  #' undocumented/secret
  K <- envelope(cells, nsim=4, saveresultof=npoints, collectrubbish=TRUE)
  #' so secret I've even forgotten how to do it
  M <- envelope(cells, nsim=4, internal=list(eject="patterns"))
}


if(ALWAYS) {
  #' Test robustness of envelope() sorting procedure when NA's are present
  #' Fails with spatstat.utils 1.12-0
  set.seed(42)
  EP <- envelope(longleaf, pcf, nsim=10, nrank=2)

  #' Test case when the maximum permitted number of failures is exceeded
  X <- amacrine[1:153] # contains exactly one point with mark='off'
  #' High probability of generating a pattern with no marks = 'off'
  E <- envelope(X, Kcross, nsim=39, maxnerr=2, maxerr.action="warn")
  A <- alltypes(X, Kcross, envelope=TRUE, nsim=39, maxnerr=2)
}

if(ALWAYS) {
  #' Internals: envelope.matrix
  Y <- matrix(rnorm(200), 10, 20)
  rr <- 1:10
  oo <- rnorm(10)
  zz <- numeric(10)
  E <- envelope(Y, rvals=rr, observed=oo, nsim=10)
  E <- envelope(Y, rvals=rr, observed=oo, jsim=1:10)
  E <- envelope(Y, rvals=rr, observed=oo, theory=zz,
                type="global", use.theory=TRUE)
  E <- envelope(Y, rvals=rr, observed=oo, theory=zz,
                type="global", use.theory=TRUE, nsim=10)
  E <- envelope(Y, rvals=rr, observed=oo, theory=zz,
                type="global", use.theory=FALSE, nsim=10)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                nsim=10, nsim2=10)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                jsim=1:10, jsim.mean=11:20)
  if(FULLTEST) print(E)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                nsim=10, jsim.mean=11:20)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                jsim=1:10, nsim2=10)
}

if(ALWAYS) {
  #' quirk with handmade summary functions ('conserve' attribute)
  Kdif <- function(X, r=NULL) { # note no ellipsis
    Y <- split(X)
    K1 <- Kest(Y[[1]], r=r)
    K2 <- Kest(Y[[2]], r=r)
    D <- eval.fv(K1-K2)
    return(D)
  }
  envelope(amacrine, Kdif, nsim=3)
}


## close 'local'
})
#
#  tests/fastK.R
#
# check fast and slow code for Kest
#       and options not tested elsewhere
#
#   $Revision: 1.5 $   $Date: 2020/04/28 12:58:26 $
#
if(ALWAYS) {
local({
  ## fast code
  Kb <- Kest(cells, nlarge=0)
  Ku <- Kest(cells, correction="none")
  Kbu <- Kest(cells, correction=c("none", "border"))
  ## slow code, full set of corrections, sqrt transformation, ratios
  Ldd <- Lest(unmark(demopat), correction="all", var.approx=TRUE, ratio=TRUE)
  ## Lotwick-Silverman var approx (rectangular window)
  Loo <- Lest(cells, correction="all", var.approx=TRUE, ratio=TRUE)
  ## Code for large dataset
  nbig <- .Machine$integer.max
  if(!is.null(nbig)) {
    nn <- ceiling(sqrt(nbig))
    if(nn < 1e6) Kbig <- Kest(runifpoint(nn),
                              correction=c("border", "bord.modif", "none"),
                              ratio=TRUE)
  }
  
  ## Kinhom
  lam <- density(cells, at="points", leaveoneout=TRUE)
  ## fast code
  Kib <- Kinhom(cells, lam, nlarge=0)
  Kiu <- Kest(cells, lam, correction="none")
  Kibu <- Kest(cells, lam, correction=c("none", "border"))
  ## slow code
  Lidd <- Linhom(unmark(demopat), sigma=bw.scott)
})

}
## 
##    tests/fvproblems.R
##
##    problems with fv, ratfv and fasp code
##
##    $Revision: 1.15 $  $Date: 2020/04/28 12:58:26 $

#' This appears in the workshop notes
#' Problem detected by Martin Bratschi

if(FULLTEST) {
local({
  Jdif <- function(X, ..., i) {
    Jidot <- Jdot(X, ..., i=i)
    J <- Jest(X, ...)
    dif <- eval.fv(Jidot - J)
    return(dif)
  }
  Z <- Jdif(amacrine, i="on")
})
}
#'
#'  Test mathlegend code
#'
local({
  K <- Kest(cells)
  if(FULLTEST) {
    plot(K)
    plot(K, . ~ r)
    plot(K, . - theo ~ r)
  }
  if(ALWAYS) {
    plot(K, sqrt(./pi)  ~ r)
  }
  if(FULLTEST) {
    plot(K, cbind(iso, theo) ~ r)
    plot(K, cbind(iso, theo) - theo ~ r)
    plot(K, sqrt(cbind(iso, theo)/pi)  ~ r)
    plot(K, cbind(iso/2, -theo) ~ r)
    plot(K, cbind(iso/2, trans/2) - theo ~ r)
  }
  if(FULLTEST) {
    ## test expansion of .x and .y
    plot(K, . ~ .x)
    plot(K, . - theo ~ .x)
    plot(K, .y - theo ~ .x)
  }
  if(ALWAYS) {
    plot(K, sqrt(.y) - sqrt(theo) ~ .x)
  }

  # problems with parsing weird strings in levels(marks(X))
  # noted by Ulf Mehlig
  if(ALWAYS) {
    levels(marks(amacrine)) <- c("Nasticreechia krorluppia", "Homo habilis")
    plot(Kcross(amacrine))
    plot(alltypes(amacrine, "K"))
  }
  if(FULLTEST) {
    plot(alltypes(amacrine, "J"))
    plot(alltypes(amacrine, pcfcross))
  }
})

#'
#'  Test quirks related to 'alim' attribute

if(FULLTEST) {
local({
  K <- Kest(cells)
  attr(K, "alim") <- NULL
  plot(K)
  attr(K, "alim") <- c(0, 0.1)
  plot(tail(K))
})
}

#'
#' Check that default 'r' vector passes the test for fine spacing

if(ALWAYS) {
local({
  a <- Fest(cells)
  A <- Fest(cells, r=a$r)
  b <- Hest(heather$coarse)
  B <- Hest(heather$coarse, r=b$r)
  # from Cenk Icos
  X <- runifpoint(100, owin(c(0,3), c(0,10)))
  FX <- Fest(X)
  FXr <- Fest(X, r=FX$r)
  JX <- Jest(X)
})
}

##' various functionality in fv.R

if(ALWAYS) {
local({
  M <- cbind(1:20, matrix(runif(100), 20, 5))
  A <- as.fv(M)
  fvlabels(A) <- c("r","%s(r)", "%s[A](r)", "%s[B](r)", "%s[C](r)", "%s[D](r)")
  A <- rename.fv(A, "M", quote(M(r)))
  A <- tweak.fv.entry(A, "V1", new.tag="r")
  A[,3] <- NULL
  A$hogwash <- runif(nrow(A))
  fvnames(A, ".") <- NULL
  #' bind.fv with qualitatively different functions
  GK <- harmonise(G=Gest(cells), K=Kest(cells))
  G <- GK$G
  K <- GK$K
  ss <- c(rep(TRUE, nrow(K)-10), rep(FALSE, 10))
  U <- bind.fv(G, K[ss, ], clip=TRUE)
  #'
  H <- rebadge.as.crossfun(K, "H", "inhom", 1, 2)
  H <- rebadge.as.dotfun(K, "H", "inhom", 3)
  #' text layout
  op <- options(width=27)
  print(K)
  options(width=18)
  print(K)
  options(op)
  #' collapse.fv
  Kb <- Kest(cells, correction="border")
  Ki <- Kest(cells, correction="isotropic")
  collapse.fv(Kb, Ki, same="theo")
  collapse.fv(anylist(B=Kb, I=Ki), same="theo")
  collapse.fv(anylist(B=Kb), I=Ki, same="theo")
  Xlist <- replicate(3, runifpoint(30), simplify=FALSE)
  Klist <- anylapply(Xlist, Kest) 
  collapse.fv(Klist, same="theo", different=c("iso", "border"))
  names(Klist) <- LETTERS[24:26]
  collapse.fv(Klist, same="theo", different=c("iso", "border"))
})
}

if(FULLTEST) {
local({
  ## rat
  K <- Kest(cells, ratio=TRUE)
  G <- Gest(cells, ratio=TRUE)
  print(K)
  compatible(K, K)
  compatible(K, G)
  H <- rat(K, attr(K, "numerator"), attr(K, "denominator"), check=TRUE)
})
}

if(FULLTEST) {
local({
  ## bug in Jmulti.R colliding with breakpts.R
  B <- owin(c(0,3), c(0,10))
  Y <- superimpose(A=runifpoint(1212, B), B=runifpoint(496, B))
  JDX <- Jdot(Y)
  JCX <- Jcross(Y)
  Jdif <- function(X, ..., i) {
    Jidot <- Jdot(X, ..., i=i)
    J <- Jest(X, ...)
    dif <- eval.fv(Jidot - J)
    return(dif)
  }
  E <- envelope(Y, Jdif, nsim=19, i="A", simulate=expression(rlabel(Y)))
})
}

if(FULLTEST) {
local({
  #' fasp axes, title, dimnames
  a <- alltypes(amacrine)
  a$title <- NULL
  plot(a, samex=TRUE, samey=TRUE)
  dimnames(a) <- lapply(dimnames(a), toupper)

  b <- as.fv(a)
})
}

if(FULLTEST) {
local({
  ## plot.anylist (fv)
  b <- anylist(A=Kcross(amacrine), B=Kest(amacrine))
  plot(b, equal.scales=TRUE, main=expression(sqrt(pi)))
  plot(b, arrange=FALSE)
})
}
