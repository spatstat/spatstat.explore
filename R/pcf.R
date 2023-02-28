#'
#'   pcf.R
#'
#'   $Revision: 1.74 $   $Date: 2023/02/28 03:14:29 $
#'
#'   calculate pair correlation function from point pattern (pcf.ppp)
#'


pcf <- function(X, ...) {
  UseMethod("pcf")
}

pcf.ppp <- function(X, ..., r=NULL,
                    kernel="epanechnikov", bw=NULL, stoyan=0.15,
                    correction=c("translate", "Ripley"),
                    divisor=c("r", "d"),
                    var.approx=FALSE,
                    domain=NULL, ratio=FALSE,
                    close=NULL)
{
  verifyclass(X, "ppp")
#  r.override <- !is.null(r)

  win <- Window(X)
  areaW <- area(win)
  npts <- npoints(X)
  samplesize <- npairs <- npts * (npts - 1)
  lambda <- npts/areaW
  lambda2area <- npairs/areaW

  kernel <- match.kernel(kernel)

  rmaxdefault <- rmax.rule("K", win, lambda)        
  
  if(!is.null(domain)) {
    # estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, win))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    # trick pcfdot() into doing it
    indom <- inside.owin(X$x, X$y, domain)
    marx <- factor(indom, levels=c(FALSE,TRUE))
    g <- pcfdot(X %mark% marx,
                i="TRUE",
                r=r,
                correction=correction, kernel=kernel, bw=bw, stoyan=stoyan,
                divisor=divisor,
                ...)
    if(!ratio) {
      ## relabel
      g <- rebadge.fv(g, quote(g(r)), "g")
    } else {
      ## construct ratfv object
      ninside <- sum(indom)
      samplesize <- ninside * (npts-1)
      g <- ratfv(as.data.frame(g), NULL, samplesize,
                 "r", quote(g(r)),
                 "theo", NULL, c(0, rmaxdefault), 
                 attr(g, "labl"), attr(g, "desc"), fname="g",
                 ratio=TRUE)
    }
    unitname(g) <- unitname(X)
    if(var.approx)
      warning("var.approx is not implemented when 'domain' is given")
    return(g)
  }

  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             good="translate",
                             best="best",
                             none="none"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, win$type, correction.given)

  divisor <- match.arg(divisor)
  
  # bandwidth
  if(is.null(bw) && (kernel == "epanechnikov")) {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(lambda)
    hmax <- h
    # conversion to standard deviation
    bw <- h/sqrt(5)
  } else if(is.numeric(bw)) {
    # standard deviation of kernel specified
    # upper bound on half-width
    hmax <- 3 * bw
  } else {
    # data-dependent bandwidth selection: guess upper bound on half-width
    hmax <- 2 * stoyan /sqrt(lambda)
  }

  ########## r values ############################
  # handle arguments r and breaks 

  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # arguments for 'density'
  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=length(r), from=0, to=rmax),
                              .StripNull = TRUE)
  
  #################################################
  
  # compute pairwise distances
  if(npts > 1) {
    needall <- any(correction %in% c("translate", "isotropic"))
    if(is.null(close)) {
      what <- if(needall) "all" else "ijd"
      close <- closepairs(X, rmax + hmax, what=what)
    } else {
      #' check 'close' has correct format
      needed <- if(!needall) c("i", "j", "d") else
                 c("i", "j", "xi", "yi", "xj", "yj", "dx", "dy", "d")
      if(any(is.na(match(needed, names(close)))))
        stop(paste("Argument", sQuote("close"),
                   "should have components named",
                   commasep(sQuote(needed))),
             call.=FALSE)
    }
    dIJ <- close$d
  } else {
    undefined <- rep(NaN, length(r))
  }

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  out <- ratfv(df,
               NULL, samplesize,
               "r", quote(g(r)),
               "theo", NULL,
               alim,
               c("r","%s[Pois](r)"),
               c("distance argument r", "theoretical Poisson %s"),
               fname="g",
               ratio=ratio)

  ###### compute #######

  bw.used <- NULL
  
  if(any(correction=="none")) {
    #' uncorrected
    if(npts > 1) {
      kdenN <- sewpcf(dIJ, 1, denargs, lambda2area, divisor)
      gN <- kdenN$g
      bw.used <- attr(kdenN, "bw")
    } else gN <- undefined
    out <- bind.ratfv(out,
                      quotient=data.frame(un=gN),
                      denominator=samplesize,
                      labl="hat(%s)[un](r)",
                      desc="uncorrected estimate of %s",
                      preferred="un",
                      ratio=ratio)
  }
  
  if(any(correction=="translate")) {
    # translation correction
    if(npts > 1) {
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=win, paired=TRUE)
      kdenT <- sewpcf(dIJ, edgewt, denargs, lambda2area, divisor)
      gT <- kdenT$g
      bw.used <- attr(kdenT, "bw")
    } else gT <- undefined
    out <- bind.ratfv(out,
                      quotient=data.frame(trans=gT),
                      denominator=samplesize,
                      labl="hat(%s)[Trans](r)",
                      desc="translation-corrected estimate of %s",
                      preferred="trans",
                      ratio=ratio)
  }

  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    if(npts > 1) {
      XI <- ppp(close$xi, close$yi, window=win, check=FALSE)
      edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
      kdenR <- sewpcf(dIJ, edgewt, denargs, lambda2area, divisor)
      gR <- kdenR$g
      bw.used <- attr(kdenR, "bw")
    } else gR <- undefined
    out <- bind.ratfv(out,
                      quotient=data.frame(iso=gR),
                      denominator=samplesize,
                      labl="hat(%s)[Ripley](r)",
                      desc="isotropic-corrected estimate of %s",
                      preferred="iso",
                      ratio=ratio)
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }

  ## variance approximation
  ## Illian et al 2008 p 234 equation 4.3.42
  if(var.approx) {
    gr <- if(any(correction == "isotropic")) gR else gT
    # integral of squared kernel
    intk2 <- kernel.squint(kernel, bw.used)
    # isotropised set covariance of window
    gWbar <- as.function(rotmean(setcov(win), result="fv"))
    vest <- gr * intk2/(pi * r * gWbar(r) * lambda^2)
    out <- bind.ratfv(out,
                      quotient=data.frame(v=vest),
                      denominator=samplesize,
                      labl="v(r)", 
                      desc="approximate variance of %s",
                      ratio=ratio)
  }

  ## Finish off
  ## default is to display all corrections
  formula(out) <- . ~ r
  fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
  ##
  unitname(out) <- unitname(X)
  ## copy to other components
  if(ratio)
    out <- conform.ratfv(out)

  attr(out, "bw") <- bw.used
  return(out)
}

# Smoothing Estimate of Weighted Pair Correlation
# d = vector of relevant distances
# w = vector of edge correction weights (in normal use)
# denargs = arguments to density.default
# lambda2area = constant lambda^2 * areaW (in normal use)

sewpcf <- function(d, w, denargs, lambda2area, divisor=c("r","d")) {
  divisor <- match.arg(divisor)
  nw <- length(w)
  if(nw != length(d) && nw != 1)
    stop("Internal error: incorrect length of weights vector in sewpcf")
  if(divisor == "d") {
    w <- w/d
    if(!all(good <- is.finite(w))) {
      nbad <- sum(!good)
      warning(paste(nbad, "infinite, NA or NaN",
                    ngettext(nbad, "contribution was", "contributions were"),
                    "deleted from pcf estimate with divisor='d'.",
                    "Fraction deleted: ",
                    paste0(round(100 * nbad/length(w), 2), "%")),
              call.=FALSE)
      d <- d[good]
      w <- w[good]
    }
    nw <- length(w)
  }
  kden <- unnormdensity(x=d, weights=w, defaults=denargs)
  r <- kden$x
  y <- kden$y 
  if(divisor == "r")
    y <- y/r
  g <- y/(2 * pi * lambda2area)
  result <- data.frame(r=r,g=g)
  attr(result, "bw") <- kden$bw
  return(result)
}

