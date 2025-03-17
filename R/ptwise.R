#'
#' ptwise.R
#'
#' Pointwise sample statistics of summary functions
#'
#'  ptwise.envelope()
#' 
#'  bias.envelope()   - pointwise bias (function of r)
#'  RMSE.envelope()   - pointwise root mean square error (function of r)
#'
#'  ISE.envelope()    - integrated squared error
#'                      (numeric vector, one value for each simulated function)
#'
#'  MISE.envelope()   - mean integrated squared error (single number)
#'  ISB.envelope()    - integrated squared bias (single number)
#'  IV.envelope()     - integrated variance (single number)
#'
#' Copyright (c) 2024-2025 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'
#' $Revision: 1.3 $ $Date: 2025/03/16 01:40:02 $
#'

ptwise.envelope <- function(object,
                            stats=c("mean", "median", "bias", 
                                    "var", "sd", "se", "mse", "rmse",
                                    "confint", "predint"),
                            ..., level=0.95, transform=NULL, theo=NULL) {
  verifyclass(object, "envelope")
  funs <- attr(object, "simfuns")
  if(is.null(funs))
    stop(paste("Envelope object does not contain",
               "the individual simulated functions;",
               "recompute the envelope with savefuns=TRUE"))
  ## apply transformation if required
  if(!is.null(transform)) {
    callenv <- parent.frame()
    funs <- eval(substitute(with(S, E, enclos=callenv),
                            list(S=quote(funs), E=transform)))
  }
  if(!is.null(theo)) stopifnot(is.function(theo))
  ## extract simulated functions as data frame
  rname <- attr(funs, "argu")
  funnames <- (names(funs) != rname)
  dd <- as.data.frame(funs)[ , funnames]
  n <- ncol(dd)
  ## handle a user-specified summary statistic
  if(is.function(stats)) {
    fx <- apply(dd, 1, stats, ...)
    values <- data.frame(r=funs$r, f=fx)
    result <- fv(values,
                 argu="r",
                 ylab=attr(funs, "ylab") %orifnull% quote(ptwise(r)),
                 valu="f",
                 fmla= f ~ r,
                 alim=attr(funs, "alim"),
                 labl=c("r", "f(r)"),
                 desc=c("distance argument r",
                        "user-supplied summary statistic"),
                 fname= attr(funs, "fname") %orifnull% "ptwise")
    attr(result, "dotnames") <- "f"
    unitname(result) <- unitname(funs)
    return(result)
  }
  
  ## compute standard summary statistics
  meanx <- medianx <- biasx <- varx <- sdx <- sex <- msex <- rmsex <- 
    lopix <- hipix <- locix <- hicix <- NULL
  stats <- match.arg(stats, several.ok=TRUE)
  if(any(hit <- (stats %in% c("bias", "mse", "rmse")))) {
    if(is.null(theo))
      stop(paste("True value 'theo' is required when stats =",
                 commasep(sQuote(stats[hit]), " or ")),
           call.=FALSE)
    truefx <- theo(funs$r)
  }
  if(any(c("mean", "confint") %in% stats))
    meanx <- apply(dd, 1, mean, na.rm=TRUE)
  if("bias" %in% stats)
    biasx <- apply(dd - truefx, 1, mean, na.rm=TRUE)
  if("median" %in% stats)
    medianx <- apply(dd, 1, median, na.rm=TRUE)
  if(any(c("var", "sd", "se", "confint") %in% stats))
    varx <- apply(dd, 1, var, na.rm=TRUE)
  if(any(c("sd", "se", "confint") %in% stats))
    sdx  <- sqrt(varx)
  if(any(c("mse", "rmse") %in% stats)) {
    msex <- apply((dd-truefx)^2, 1, mean, na.rm=TRUE)
    if("rmse" %in% stats)
      rmsex <- sqrt(msex)
  }
  if(any(c("se", "confint") %in% stats))
    sex  <- sdx/sqrt(n)
  if(any(c("confint", "predint") %in% stats)) {
    if(level > 1) stop("level should be a fraction between 0 and 1")
    ahi <- (1+level)/2
    alo <- (1-level)/2
    if("confint" %in% stats) {
      df <- ncol(dd) - 1
      locix <- meanx + qt(ahi, df) * sex
      hicix <- meanx + qt(alo, df) * sex
    }
    if("predint" %in% stats) {
      lopix <- apply(dd, 1, quantile, probs=alo, na.rm=TRUE)
      hipix <- apply(dd, 1, quantile, probs=ahi, na.rm=TRUE)
    }
  }
  ## now remove statistics that were not requested
  if(!("mean" %in% stats)) meanx <- NULL 
  if(!("var" %in% stats)) varx <- NULL
  if(!("sd" %in% stats)) sdx <- NULL
  if(!("se" %in% stats)) sex <- NULL
  if(!("mse" %in% stats)) msex <- NULL 
 
  ## assemble results
  values <- as.data.frame(cbind(
    r      = funs$r,
    mean   = meanx,
    median = medianx,
    bias   = biasx,
    var    = varx,
    sd     = sdx,
    se     = sex,
    mse    = msex,
    rmse   = rmsex,
    loci   = locix,
    hici   = hicix,
    lopi   = lopix,
    hipi   = hipix))

  Labels <- c(
    r      = "r",
    mean   = "mean(r)",
    median = "median(r)",
    bias   = "bias(r)",
    var    = "var(r)",
    sd     = "sd(r)",
    mse    = "mse(r)",
    rmse   = "rmse(r)",
    loci   = "lo[CI](r)",
    hici   = "hi[CI](r)",
    lopi   = "lo[PI](r)",
    hipi   = "hi[PI](r)")

  Descrips <- c(
    r      = "distance argument r",
    mean   = "pointwise mean",
    median = "pointwise median",
    bias   = "pointwise bias",
    var    = "pointwise variance",
    sd     = "pointwise standard deviation",
    se     = "pointwise standard error",
    mse    = "pointwise mean squared error",
    rmse   = "pointwise root-mean-squared error",
    loci   = "lower limit of confidence interval for true mean",
    hici   = "upper limit of confidence interval for true mean",
    lopi   = "lower limit of prediction interval for function value",
    hipi   = "upper limit of prediction interval for function value")

  chosen <- colnames(values)
  labl <- unname(Labels[chosen])
  desc <- unname(Descrips[chosen])

  statlabels <- setdiff(chosen, "r")
  
  result <- fv(values,
               argu="r",
               ylab=attr(funs, "ylab") %orifnull% quote(ptwise(r)),
               fmla= . ~ r,
               valu=statlabels[1],
               alim=attr(funs, "alim"),
               labl=labl,
               desc=desc,
               fname = attr(funs, "fname") %orifnull% "ptwise")
  attr(result, "dotnames") <- statlabels

  if(any(c("confint", "predint") %in% stats)) {
    fvnames(result, ".s") <-
      if("confint" %in% stats) c("loci", "hici") else c("lopi", "hipi")
  }
  
  unitname(result) <- unitname(funs)
  return(result)
}


bias.envelope <- function(object, theo, CI=TRUE, level=0.95) {
  if(!CI) {
    a <- ptwise.envelope(object, "bias", theo=theo)
    return(a)
  }
  a <- ptwise.envelope(object, c("mean", "confint"))
  stopifnot(is.function(theo))
  a[] <- with(a, . - theo(.x))[drop=TRUE]
  colnames(a) <- sub("mean", "bias", colnames(a))
  attr(a, "labl") <- sub("mean", "bias", attr(a, "labl"))
  attr(a, "desc") <- sub("mean", "bias", attr(a, "desc"))
  a <- rebadge.fv(a, new.fname="bias", new.ylab=quote(bias(r)))
  return(a)
}

RMSE.envelope <- function(object, theo) {
  f <- ptwise.envelope(object, "rmse", theo=theo)
  rebadge.fv(f, new.yexp=quote(RMSE(r)), new.fname="RMSE")
}

#'  Integrated squared error
#'  - one numerical value for each simulated function

ISE.envelope <- function(object, theo, domain=NULL, dimension=2) {
  ## Integrated squared error
  stopifnot(is.function(theo))
  verifyclass(object, "envelope")
  ## extract function estimates
  funs <- attr(object, "simfuns")
  if(is.null(funs))
    stop(paste("Envelope object does not contain",
               "the individual simulated functions;",
               "recompute the envelope with savefuns=TRUE"))
  ## trim to domain
  if(!is.null(domain)) {
    check.range(domain)
    rr <- with(funs, domain[1] <= .x & domain[2] >= .x)
    funs <- funs[rr, , drop=FALSE]
  }
  ## extract function values
  df <- as.data.frame(funs)
  i <- match(fvnames(funs, ".x"), colnames(df))
  fx <- df[ , -i]
  x  <- df[ , i, drop=TRUE]
  ## squared error
  y <- (fx - theo(x))^2
  ## integrate
  if(dimension == 2)
    y <- 2 * pi * x * y
  n <- nrow(y)
  z <- colSums(diff(x) * (y[ -1, ,drop=FALSE] + y[-nrow(y), ,drop=FALSE]))/2
  return(z)
}

#' Integrated moments - single numerical value returned

MISE.envelope <- function(object, theo, domain, dimension=2) {
  ## mean integrated squared error
  mse <- ptwise.envelope(object, "mse", theo=theo)
  stopifnot(dimension %in% c(1,2))
  if(dimension == 2)
    mse <- with(mse, 2 * pi * .x * .a)
  unname(integral(mse, domain=domain))
}

ISB.envelope <- function(object, theo, domain, dimension=2) {
  ## integrated squared bias
  b <- ptwise.envelope(object, "bias", theo=theo)
  bs <- eval.fv(b^2)
  stopifnot(dimension %in% c(1,2))
  if(dimension == 2)
    bs <- with(bs, 2 * pi * .x * .a)
  unname(integral(bs, domain=domain))
}

IV.envelope <- function(object, domain, dimension=2) {
  ## integrated variance
  b <- ptwise.envelope(object, "var")
  stopifnot(dimension %in% c(1,2))
  if(dimension == 2)
    b <- with(b, 2 * pi * .x * .a)
  unname(integral(b, domain=domain))
}

