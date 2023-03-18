##
##   Math.fv.R
##
##   Inline arithmetic for 'fasp'
##
##   $Revision: 1.1 $ $Date: 2023/03/18 07:06:11 $


Math.fasp <- function(x, ...){
  eval(substitute(eval.fasp(G(x)),
                  list(G=as.name(.Generic),
                       x=substitute(x))))
}

Complex.fasp <- function(z){
  eval(substitute(eval.fasp(G(z)),
                  list(G=as.name(.Generic),
                       z=substitute(z))))
}

Ops.fasp <- function(e1,e2=NULL) {
  eval(substitute(eval.fasp(G(e1,e2)), list(G=as.name(.Generic),
                                            e1=substitute(e1),
                                            e2=substitute(e2))))
}

Summary.fasp <- local({
  
  Summary.fasp <- function(..., na.rm=FALSE){
    argh <- list(...)
    arrays <- sapply(argh, inherits, what="fasp")
    argh[arrays] <- lapply(argh[arrays], processArray, op=.Generic, na.rm=na.rm)
    funs <- sapply(argh, is.fv)
    if(any(funs)) 
      argh[funs] <- lapply(argh[funs], .Generic, na.rm=na.rm)
    do.call(.Generic, c(argh, list(na.rm = na.rm)))
  }

  processArray <- function(x, op, na.rm=FALSE) {
    ## extract individual fv objects and apply operation 'op'
    y <- unlist(lapply(x$fns, op, na.rm=na.rm))
    ## apply 'op' to the results
    do.call(op, c(y, list(na.rm=na.rm)))
  }
      
  Summary.fasp
})

