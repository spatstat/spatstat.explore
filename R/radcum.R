##
## radcum.R
##
## cumulative integral as function of distance
##
##  $Revision: 1.1 $ $Date: 2022/07/16 03:29:09 $

radcumint <- function(X, ..., origin, Xname, result=c("fv", "im")) {
  if(missing(Xname))
    Xname <- sensiblevarname(short.deparse(substitute(X)), "X")
  result <- match.arg(result)
  trap.extra.arguments(..., .Context="radcum")
  stopifnot(is.im(X))
  if(!missing(origin) && !is.null(origin)) {
    X <- shift(X, origin=origin)
    backshift <- -getlastshift(X)
  } else {
    backshift <- NULL
  }
  #' determine discretisation steps
  rmax <- with(vertices(Frame(X)), sqrt(max(x^2+y^2)))
  pixarea <- with(X, xstep * ystep)
  eps <- with(X, sqrt(xstep^2 + ystep^2))
  #'
  Xdata <- as.data.frame(X)
  values <- Xdata$value
  radii <- with(Xdata, sqrt(x^2+y^2))
  #'
  rmax <- max(max(radii), rmax) + 2*eps
  rr <- seq(0, rmax, by=eps)
  wh <- whist(radii, breaks=rr, weights=pixarea * values)
  yy <- cumsum(c(0,wh))
  #'
  df <- data.frame(r=rr, f=yy)
  FUN <- fv(df,
            argu="r",
            ylab=substitute(bar(X)(r), list(X=as.name(Xname))),
            valu="f",
            fmla=(. ~ r),
            alim=c(0, rmax),
            labl=c("r", "%s(r)"),
            desc=c("distance argument r",
                   "rotational average"),
            unitname=unitname(X),
            fname=paste0("bar", paren(Xname)))
  fvnames(FUN, ".") <- "f"
  if(result == "fv") return(FUN)
  ## compute image
  FUN <- as.function(FUN)
  IM <- as.im(function(x,y,FUN){ FUN(sqrt(x^2+y^2)) }, X, FUN=FUN)
  if(!is.null(backshift))
    IM <- shift(IM,backshift)
  return(IM)
}
