##
##   Math.fv.R
##
##   Inline arithmetic for 'fv' 
##
##   $Revision: 1.7 $ $Date: 2023/03/18 10:14:26 $


Math.fv <- function(x, ...){
  force(x)
  eval(substitute(eval.fv(G(x)),
                  list(G=as.name(.Generic),
                       x=quote(x))))
}

Complex.fv <- function(z){
  force(z)
  eval(substitute(eval.fv(G(z)),
                  list(G=as.name(.Generic),
                       z=quote(z))))
}

Ops.fv <- function(e1,e2=NULL) {
  m <- match.call()
  if(is.name(m$e1) || (is.atomic(m$e1) && length(m$e1) == 1)) {
    e1use <- substitute(e1)
  } else {
    force(e1)
    e1use <- quote(e1)
  }
  if(is.name(m$e2) || (is.atomic(m$e2) && length(m$e2) == 1)) {
    e2use <- substitute(e2)
  } else {
    force(e2)
    e2use <- quote(e2)
  }
  eval(substitute(eval.fv(G(e1,e2)), list(G=as.name(.Generic),
                                          e1=e1use,
                                          e2=e2use)))
}

Summary.fv <- local({
  
  Summary.fv <- function(..., na.rm=FALSE){
    argh <- list(...)
    funs <- sapply(argh, is.fv)
    argh[funs] <- lapply(argh[funs], getValues)
    do.call(.Generic, c(argh, list(na.rm = na.rm)))
  }

  getValues <- function(x) {
    xdat <- as.matrix(as.data.frame(x))
    yall <- fvnames(x, ".")
    vals <- xdat[, yall]
    return(as.vector(vals))
  }
  
  Summary.fv
})


