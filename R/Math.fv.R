##
##   Math.fv.R
##
##   Inline arithmetic for 'fv' 
##
##   $Revision: 1.5 $ $Date: 2023/03/18 07:06:44 $


Math.fv <- function(x, ...){
  eval(substitute(eval.fv(G(x)),
                  list(G=as.name(.Generic),
                       x=substitute(x))))
}

Complex.fv <- function(z){
  eval(substitute(eval.fv(G(z)),
                  list(G=as.name(.Generic),
                       z=substitute(z))))
}

Ops.fv <- function(e1,e2=NULL) {
  eval(substitute(eval.fv(G(e1,e2)), list(G=as.name(.Generic),
                                          e1=substitute(e1),
                                          e2=substitute(e2))))
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


