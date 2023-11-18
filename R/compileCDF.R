#'     compileCDF.R
#'
#'     Wrappers for estimating CDF from right-censored data
#'
#'   Copyright (c) 1991-2023 Adrian Baddeley, Rolf Turner and Ege Rubak
#'   GNU Public Licence (>= 2.0)
#' 
#'   $Revision: 1.1 $ $Date: 2023/11/04 04:45:32 $

censtimeCDFest <- function(o, cc, d, breaks, ...,
                           KM=TRUE, RS=TRUE, HAN=TRUE, RAW=TRUE,
                           han.denom=NULL, tt=NULL, pmax=0.9,
                           fname="CDF", fexpr=quote(CDF(r))) {
# Histogram-based estimation of cumulative distribution function
# of lifetimes subject to censoring.
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#       han.denom: denominator (eroded area) for each value of r
#       tt: uncensored lifetimes T_i, if known  
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
  rval <- breaks$r
  rmax <- breaks$max
  # Kaplan-Meier and/or Reduced Sample
  out <- km.rs.opt(o, cc, d, breaks, KM=KM, RS=RS)
  # convert to data frame
  out$breaks <- NULL
  df <- as.data.frame(out)
  # Raw ecdf of observed lifetimes if available
  if(RAW && !is.null(tt)) {
    h <- whist(tt[tt <= rmax], breaks=bval)
    df <- cbind(df, data.frame(raw=cumsum(h)/length(tt)))
  }
  # Hanisch
  if(HAN) {
    if(is.null(han.denom))
      stop("Internal error: missing denominator for Hanisch estimator")
    if(length(han.denom) != length(rval))
      stop(paste("Internal error:",
                 "length(han.denom) =", length(han.denom),
                 "!=", length(rval), "= length(rvals)"))
    #  uncensored distances
    x <- o[d]
    # calculate Hanisch estimator
    h <- whist(x[x <= rmax], breaks=bval)
    H <- cumsum(h/han.denom)
    df <- cbind(df, data.frame(han=H/max(H[is.finite(H)])))
  }
  # determine appropriate plotting range
  bestest <- if(KM) "km" else if(HAN) "han" else if(RS) "rs" else "raw"
  alim <- range(df$r[df[[bestest]] <= pmax])
  # convert to fv object
  nama <-  c("r",  "km", "hazard", "han", "rs", "raw")
  avail <- c(TRUE,  KM,  KM,       HAN,   RS,   RAW)
  iscdf <- c(FALSE, TRUE, FALSE,   TRUE,  TRUE, TRUE)
  labl <- c("r",
            makefvlabel(NULL, "hat", fname, "km"),
            "hat(lambda)(r)",
            makefvlabel(NULL, "hat", fname, "han"),
            makefvlabel(NULL, "hat", fname, "bord"),
            makefvlabel(NULL, "hat", fname, "raw")
            )[avail]
  desc <- c("distance argument r",
            "Kaplan-Meier estimate of %s",
            "Kaplan-Meier estimate of hazard function lambda(r)",
            "Hanisch estimate of %s",
            "border corrected estimate of %s",
            "uncorrected estimate of %s")[avail]
  df <- df[, nama[avail]]
  Z <- fv(df, "r", fexpr, bestest, . ~ r, alim, labl, desc,
          fname=fname)
  fvnames(Z, ".") <- nama[iscdf & avail]
  return(Z)
}

# simple interface for students and code development

compileCDF <- function(D, B, r, ..., han.denom=NULL, check=TRUE) {
  han <- !is.null(han.denom)
  breaks <- breakpts.from.r(r)
  if(check) {
    stopifnot(length(D) == length(B) && all(D >= 0) && all(B >= 0))
    if(han)
      stopifnot(length(han.denom) == length(r))
  }
  D <- as.vector(D)
  B <- as.vector(B)
  # observed (censored) lifetimes
  o <- pmin.int(D, B)
  # censoring indicators
  d <- (D <= B)
  # go
  result <- censtimeCDFest(o, B, d, breaks,
                           HAN=han, 
                           han.denom=han.denom,
                           RAW=TRUE, tt=D)
  result <- rebadge.fv(result, new.fname="compileCDF")
}
