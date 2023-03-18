##
##   Math.fv.R
##
##   Inline arithmetic for 'fasp'
##
##   $Revision: 1.3 $ $Date: 2023/03/18 10:14:35 $


Math.fasp <- function(x, ...){
  force(x)
  eval(substitute(eval.fasp(G(x)),
                  list(G=as.name(.Generic),
                       x=quote(x))))
}

Complex.fasp <- function(z){
  force(z)
  eval(substitute(eval.fasp(G(z)),
                  list(G=as.name(.Generic),
                       z=quote(z))))
}

Ops.fasp <- function(e1,e2=NULL) {
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
  eval(substitute(eval.fasp(G(e1,e2)), list(G=as.name(.Generic),
                                            e1=e1use,
                                            e2=e2use)))
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

