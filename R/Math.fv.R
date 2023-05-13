##
##   Math.fv.R
##
##   Inline arithmetic for 'fv' 
##
##   $Revision: 1.9 $ $Date: 2023/05/13 01:11:16 $


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
  objects <- list()
  if(is.name(m$e1) || (is.atomic(m$e1) && length(m$e1) == 1)) {
    ## e1 is the name of an fv object, or is a single value.
    ## It will appear directly in the resulting function name
    e1use <- substitute(e1)
  } else {
    ## e1 is an expression that should first be evaluated
    ## It will appear as 'e1' in the resulting function name
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
  eval(substitute(eval.fv(G(e1,e2),
                          envir=evalframe),
                  list(G=as.name(.Generic),
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


