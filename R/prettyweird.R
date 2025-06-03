#'
#' prettyweird.R
#'
#' 
#' prettyweird(f) where f is a CDF
#' or
#' prettyweird(x, p) where x is a vector of quantiles
#'                   and p is the vector of corresponding probabilities.
#'

prettyweird <- function(x, p, pieces=3, neach=2) {
  #' find a 'pretty' sequence of QUANTILES to be plotted on probability scale.
  if(missing(p) && inherits(x, c("ecdf", "ewcdf"))) {
    ## x is a CDF
    Fx <- x
    Finv <- quantilefun(x)
  } else {
    check.nvector(x)
    check.nvector(p)
    stopifnot(length(x) == length(p))
    Fx <- approxfun(x, p, yleft=0, yright=1)
    Finv <- approxfun(p, x, rule=2)
  }
  breaks <- Finv((0:pieces)/pieces)
  XX <- numeric(0)
  for(i in seq_len(pieces)) {
    XX <- c(XX, pretty(breaks[c(i, i+1)], n=neach+2, min.n=neach))
  }
  XX <- sort(unique(XX))
  if(any(close <- diff(XX) < 0.02 * mean(diff(XX))))
    XX <- XX[c(TRUE, !close)]
  PP <- Fx(XX)
  return(data.frame(x=XX, p=PP))
}
