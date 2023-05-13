##
##   Math.fv.R
##
##   Inline arithmetic for 'fasp'
##
##   $Revision: 1.4 $ $Date: 2023/05/13 01:11:06 $


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
  objects <- list()
  if(is.name(m$e1) || (is.atomic(m$e1) && length(m$e1) == 1)) {
    ## e1 is the name of a fasp object, or is a single value.
    ## It will appear directly in the resulting function names
    e1use <- substitute(e1)
  } else {
    ## e1 is an expression that should first be evaluated
    ## It will appear as 'e1' in the resulting function names
    e1use <- quote(e1)
    objects$e1 <- eval(e1)
  }
  if(is.name(m$e2) || (is.atomic(m$e2) && length(m$e2) == 1)) {
    e2use <- substitute(e2)
  } else {
    e2use <- quote(e2)
    objects$e2 <- eval(e2)
  }
  callframe <- parent.frame()
  evalframe <-
    if(length(objects)) list2env(objects, parent=callframe) else callframe
  eval(substitute(eval.fasp(G(e1,e2),
                            envir=evalframe),
                  list(G=as.name(.Generic),
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

