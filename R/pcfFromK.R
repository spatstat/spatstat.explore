#'
#' pcfFromK.R
#' 
#' Calculate pcf from other estimators of K or Kcross 
#'
#'     pcf.fv, pcf.fasp
#'
#' $Revision: 1.1 $ $Date: 2023/02/19 01:02:57 $


pcf.fasp <- function(X, ..., method="c") {
  verifyclass(X, "fasp")
  Y <- X
  Y$title <- paste("Array of pair correlation functions",
                   if(!is.null(X$dataname)) "for",
                   X$dataname)
  # go to work on each function
  for(i in seq_along(X$fns)) {
    Xi <- X$fns[[i]]
    PCFi <- pcf.fv(Xi, ..., method=method)
    Y$fns[[i]] <- PCFi
    if(is.fv(PCFi))
      Y$default.formula[[i]] <- formula(PCFi)
  }
  return(Y)
}


pcf.fv <- local({

  callmatched <- function(fun, argue) {
    formalnames <- names(formals(fun))
    formalnames <- formalnames[formalnames != "..."]
    do.call(fun, argue[names(argue) %in% formalnames])
  }

  pcf.fv <- function(X, ..., method="c") {
    verifyclass(X, "fv")
  
    # extract r and the recommended estimate of K
    r <- with(X, .x)
    K <- with(X, .y)
    alim <- attr(X, "alim")

    # remove NA's
    ok <- !is.na(K)
    K <- K[ok]
    r <- r[ok]
    switch(method,
           a = {
             ss <- callmatched(smooth.spline,
                               list(x=r, y=K, ...))
             dK <- predict(ss, r, deriv=1)$y
             g <- dK/(2 * pi * r)
           },
           b = {
             y <- K/(2 * pi * r)
             y[!is.finite(y)] <- 0
             ss <- callmatched(smooth.spline,
                               list(x=r, y=y, ...))
             dy <- predict(ss, r, deriv=1)$y
             g <- dy + y/r
           },
           c = {
             z <- K/(pi * r^2)
             z[!is.finite(z)] <- 1
             ss <- callmatched(smooth.spline,
                               list(x=r, y=z, ...))
             dz <- predict(ss, r, deriv=1)$y
             g <- (r/2) * dz + z
           },
           d = {
             z <- sqrt(K)
             z[!is.finite(z)] <- 0
             ss <- callmatched(smooth.spline,
                               list(x=r, y=z, ...))
             dz <- predict(ss, r, deriv=1)$y
             g <- z * dz/(pi * r)
           },
           stop(paste("unrecognised method", sQuote(method)))
           )

    # pack result into "fv" data frame
    Z <- fv(data.frame(r=r,
                       theo=rep.int(1, length(r)),
                       pcf=g),
            "r", substitute(g(r), NULL), "pcf", . ~ r, alim,
            c("r", "%s[pois](r)", "%s(r)"),
            c("distance argument r",
              "theoretical Poisson value of %s",
              "estimate of %s by numerical differentiation"),
            fname="g")
    unitname(Z) <- unitname(X)
    return(Z)
  }

  pcf.fv
})

