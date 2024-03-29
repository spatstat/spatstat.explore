#
#     eval.fasp.R
#
#
#        eval.fasp()             Evaluate expressions involving fasp objects
#
#        compatible.fasp()       Check whether two fasp objects are compatible
#
#     $Revision: 1.13 $     $Date: 2023/03/18 10:25:19 $
#

eval.fasp <- local({

  eval.fasp <- function(expr, envir, dotonly=TRUE) {
    #' convert syntactic expression to 'expression' object
    e <- as.expression(substitute(expr))
    #' convert syntactic expression to call
    ##  elang <- substitute(expr)
    #' find names of all variables in the expression
    varnames <- all.vars(e)
    if(length(varnames) == 0)
      stop("No variables in this expression")
    ## get the actual variables
    if(missing(envir)) {
      envir <- sys.parent()
    } else if(is.list(envir)) {
      envir <- list2env(envir, parent=parent.frame())
    }
    vars <- lapply(as.list(varnames), get, envir=envir)
    names(vars) <- varnames
    ## find out which ones are fasp objects
    isfasp <- unlist(lapply(vars, inherits, what="fasp"))
    if(!any(isfasp))
      stop("No fasp objects in this expression")
    fasps <- vars[isfasp]
    nfasps <- length(fasps)
    ## test whether the fasp objects are compatible
    if(nfasps > 1L && !(do.call(compatible, unname(fasps))))
      stop(paste(if(nfasps > 2) "some of" else NULL,
                 "the objects",
                 commasep(sQuote(names(fasps))),
                 "are not compatible"))
    ## copy first object as template
    result <- fasps[[1L]]
    which <- result$which
    nr <- nrow(which)
    nc <- ncol(which)
    ## create environment for evaluation
    fenv <- new.env()
    ## for each [i,j] extract fv objects and evaluate expression
    for(i in seq_len(nr))
      for(j in seq_len(nc)) {
        ## extract fv objects at position [i,j]
        funs <- lapply(fasps, getpanel, i=i, j=j)
        ## insert into list of argument values
        vars[isfasp] <- funs
        ## assign them into the right environment
        for(k in seq_along(vars)) 
          assign(varnames[k], vars[[k]], envir=fenv)
        ## evaluate
        resultij <- eval(substitute(eval.fv(ee,ff,dd),
                                    list(ee=e[[1]], ff=fenv, dd=dotonly)))
        ## insert back into fasp
        result$fns[[which[i,j] ]] <- resultij
      }
    result$title <- paste("Result of eval.fasp(", e, ")", sep="")
    return(result)
  }

  getpanel <- function(x, i, j) { as.fv(x[i,j]) }

  eval.fasp
})

compatible.fasp <- function(A, B, ...) {
  verifyclass(A, "fasp")
  if(missing(B)) return(TRUE)
  verifyclass(B, "fasp")
  dimA <- dim(A$which)
  dimB <- dim(B$which)
  if(!all(dimA == dimB))
    return(FALSE)
  for(i in seq_len(dimA[1L])) 
    for(j in seq_len(dimA[2L])) {
      Aij <- as.fv(A[i,j])
      Bij <- as.fv(B[i,j])
      if(!compatible.fv(Aij, Bij))
        return(FALSE)
    }
  # A and B agree
  if(length(list(...)) == 0) return(TRUE)
  # recursion
  return(compatible.fasp(B, ...))
}

