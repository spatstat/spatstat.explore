# compileK
#
# Function to take a matrix of pairwise distances
# and compile a 'K' function in the format required by spatstat.
#
#   $Revision: 1.16 $  $Date: 2023/08/16 02:07:59 $
# -------------------------------------------------------------------

compileK <- function(D, r, weights=NULL, denom=1, check=TRUE, ratio=FALSE,
                     fname="K", samplesize=denom) {
  # process r values
  breaks <- breakpts.from.r(r)
  rmax <- breaks$max
  r    <- breaks$r
  # check that D is a symmetric matrix with nonnegative entries
  if(check)
    stopifnot(is.matrix(D) && isSymmetric(D) && all(D >= 0))
  # ignore the diagonal; throw away any D values greater than rmax
  ok <- (D <= rmax & D > 0)
  Dvalues <- D[ok]
  #
  # weights?
  if(!is.null(weights)) {
    stopifnot(is.matrix(weights) && all(dim(weights)==dim(D)))
    wvalues <- weights[ok]
  } else wvalues <- NULL
  # count the number of D values in each interval (r[k], r[k+1L]]
  counts <- whist(Dvalues, breaks=breaks$val, weights=wvalues)
  # cumulative counts: number of D values in [0, r[k])
  Kcount <- cumsum(counts)
  # calculate estimate
  Kratio <- Kcount/denom
  # wrap it up as an 'fv' object for use in spatstat
  df <- data.frame(r=r, est=Kratio)
  labl <- c("r", makefvlabel(NULL, "hat", fname))
  K <- fv(df, "r", quote(K(r)), "est", . ~ r , c(0,rmax), labl,
          c("distance argument r", "estimated %s"),
          fname=fname)
  if(ratio) {
    if(missing(samplesize) || is.null(samplesize)) {
      Numer <- Kcount
      Denom <- denom
    } else {
      ## adjust numer/denom so that denominator is sample size
      Numer <- Kcount * samplesize/denom
      Denom <- samplesize
    }
    ## create numerator and denominator as fv objects
    Knum <- fv(data.frame(r=r, est=Numer),
               "r", quote(K(r)), "est", . ~ r , c(0,rmax), labl,
               c("distance argument r", "numerator of estimated %s"),
               fname=fname)
    Kden <- fv(data.frame(r=r, est=Denom),
               "r", quote(K(r)), "est", . ~ r , c(0,rmax), labl,
               c("distance argument r", "denominator of estimated %s"),
               fname=fname)
    K <- rat(K, Knum, Kden, check=FALSE)
  }
  return(K)
}


compilepcf <- function(D, r, weights=NULL, denom=1, check=TRUE,
                       endcorrect=TRUE, ratio=FALSE, ...,
                       fname="g", samplesize=denom) {
  # process r values
  breaks <- breakpts.from.r(r)
  if(!breaks$even)
    stop("compilepcf: r values must be evenly spaced", call.=FALSE)
  r    <- breaks$r
  rmax <- breaks$max
  # check that D is a symmetric matrix with nonnegative entries
  if(check)
    stopifnot(is.matrix(D) && isSymmetric(D) && all(D >= 0))
  # ignore the diagonal; throw away any D values greater than rmax
  ok <- (D <= rmax & D > 0)
  Dvalues <- D[ok]
  #
  # weights?
  if(!is.null(weights)) {
    stopifnot(is.matrix(weights) && all(dim(weights)==dim(D)))
    wvalues <- weights[ok]
    totwt <- sum(wvalues)
    normwvalues <- wvalues/totwt
  } else {
    nv <- length(Dvalues)
    normwvalues <- rep.int(1/nv, nv)
    totwt <- nv
  }
  # form kernel estimate
  rmin <- min(r)
  rmax <- max(r)
  nr   <- length(r)
  Ddens <- do.call.matched(density.default,
                           resolve.defaults(
                             list(x=Dvalues,
                                  weights=normwvalues,
                                  from=rmin, to=rmax, n=nr),
                             list(...),
                             list(warnWbw=FALSE)))
  gval <- Ddens$y * totwt
  # normalise
  gval <- gval/denom
  # edge effect correction at r = 0
  if(endcorrect) {
    one <- do.call.matched(density.default,
                           resolve.defaults(
                             list(x=seq(rmin,rmax,length=512),
                                  bw=Ddens$bw,
                                  adjust=1,
                                  from=rmin, to=rmax, n=nr),
                             list(...),
                             list(warnWbw=FALSE)
                           ))
    onefun <- approxfun(one$x, one$y, rule=2)
    gval <- gval /((rmax-rmin) * onefun(Ddens$x))
  }
  # wrap it up as an 'fv' object for use in spatstat
  df <- data.frame(r=r, est=gval)
  if(!ratio) {
    g <- fv(df, "r", quote(g(r)), "est", . ~ r , c(0,rmax),
            c("r", makefvlabel(NULL, "hat", fname)),
    	    c("distance argument r", "estimated %s"),
	    fname=fname)
  } else {
    if(is.null(samplesize)) samplesize <- denom
    num <- data.frame(r=r, est=gval * samplesize)
    den <- data.frame(r=r, est=samplesize)
    g <- ratfv(df=NULL, numer=num, denom=den,
               "r", quote(g(r)), "est", . ~ r , c(0,rmax),
               c("r", makefvlabel(NULL, "hat", fname)), 
               c("distance argument r", "estimated %s"),
               fname=fname)
  }
  attr(g, "bw") <- Ddens$bw
  return(g)
}
