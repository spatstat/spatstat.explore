#' nncount.R
#'
#' Functions nncount() and nnequal()
#' 
#' Copyright (c) 2019-2025 Lucia Cobo Sanchez and Adrian Baddeley
#' GNU Public Licence >= 2.0
#'
#' $Revision: 1.1 $ $Date: 2025/09/28 03:11:58 $

nncount <- function(X, i=1, j=2, ..., kmax=20, ratio=TRUE, cumulative=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  marx <- marks(X)
  lev <- levels(marx)
  if(is.numeric(i)) i <- lev[i]
  if(is.numeric(j)) j <- lev[j]
  if(is.na(match(i, lev)))
    stop(paste("Unrecognised value", i, "for argument i"))
  if(is.na(match(j, lev)))
    stop(paste("Unrecognised value", j, "for argument j"))
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  N <- unname(nnwhich(X, k=1:kmax))
  mi <- (marx == i)
  mj <- (marx == j)
  if(!any(mi)) {
    numer <- denom <- rep(0, kmax)
  } else if(!any(mj)) {
    numer <- rep(0, kmax)
    denom <- sum(mi)
  } else {
    Xsub <- X[mi]
    Nsub <- N[mi, , drop=FALSE]
    Dsub <- unname(nndist(Xsub,  k=1:kmax))
    Bsub <- bdist.points(Xsub)
    observed <- (Dsub <= Bsub)
    counted <- matrix(mj[Nsub], ncol=kmax)
    if(cumulative) {
      numer <- rowSums(apply(observed & counted, 1, cumsum))
      denom <- rowSums(apply(observed,           1, cumsum))
    } else {
      numer <- colSums(observed & counted)
      denom <- colSums(observed)
    }
  }
  estimate <- ifelse(denom > 0, numer/denom, 0)

  pj <- mean(mj)
  df <- data.frame(k=1:kmax, theo=pj, bord=estimate)
  desc <- c("neighbour order k",
            "theoretical %s", 
            "border-corrected estimate of %s")
  labl <- c("k","{%s[%s]^{theo}}(k)", "{hat(%s)[%s]^{bord}}(k)")
  dendf <- data.frame(k=1:kmax, theo=denom, bord=denom)
  ylab <- substitute(N[i,j](r), list(i=iname,j=jname))
  yexp <- substitute(N[list(i,j)](k), list(i=iname,j=jname))
  fname <- c("N", paste0("list(", iname, ",", jname, ")"))
  Z <- ratfv(df, NULL, dendf, 
             argu="k",
             ylab=ylab,
             valu="bord",
             fmla = . ~ k,
             alim=c(1,kmax),
             labl=labl, desc=desc, fname=fname, yexp=yexp,
             ratio=ratio,
             unitname=c("neighbour step", "neighbour steps"))
  ##
  return(Z)
}

nnequal <- function(X, ..., kmax=20, ratio=TRUE, cumulative=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  N <- unname(nnwhich(X, k=1:kmax))
  D <- unname(nndist(X,  k=1:kmax))
  B <- bdist.points(X)
  observed <- (D <= B)
  marx <- marks(X)
  mI <- matrix(marx[row(N)], ncol=kmax)
  mJ <- matrix(marx[N], ncol=kmax)
  counted <- (mI == mJ)
  if(cumulative) {
    numer <- rowSums(apply(observed & counted, 1, cumsum))
    denom <- rowSums(apply(observed,           1, cumsum))
  } else {
    numer <- colSums(observed & counted)
    denom <- colSums(observed)
  }
  estimate <- ifelse(denom > 0, numer/denom, 0)

  m <- as.integer(table(marx))
  n <- npoints(X)
  pequal <- sum(m*(m-1)/(n*(n-1)))
  
  df <- data.frame(k=1:kmax, theo=pequal, bord=estimate)
  desc <- c("neighbour order k",
            "theoretical %s", 
            "border-corrected estimate of %s")
  labl <- c("k","%s[theo](k)", "hat(%s)[bord](k)")
  dendf <- data.frame(k=1:kmax, theo=denom, bord=denom)
  yexp <- ylab <- quote(E(k))
  fname <- "E"
  Z <- ratfv(df, NULL, dendf, 
             argu="k",
             ylab=ylab,
             valu="bord",
             fmla = . ~ k,
             alim=c(1,kmax),
             labl=labl, desc=desc, fname=fname, yexp=yexp,
             ratio=ratio,
             unitname=c("neighbour step", "neighbour steps"))
  ##
  return(Z)
}

