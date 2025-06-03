#'  Youden statistic
#'  (max deviation from diagonal)

youden <- function(X, sign=c("positive", "absolute", "negative")) {
  verifyclass(X, "roc")
  sign <- match.arg(sign)
  ## choose meaningful columns
  wanted <- setdiff(fvnames(X, "."), "null")
  ## compute difference from null
  df <- as.data.frame(X)
  differ <- df[, wanted, drop=FALSE] - df[,"null"]
  ## calculate max deviation
  deviation <- switch(sign,
                      positive = differ,
                      absolute = abs(differ),
                      negative = -differ)
  y <- apply(deviation, 2, max)
  return(pmax(y,0))
}
