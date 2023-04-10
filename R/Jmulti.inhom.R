#'
#'    Jmulti.inhom.R
#'
#'    Inhomogeneous multitype J function
#'
#'    original code by Jonatan Gonzalez
#'    Edited for spatstat by Adrian Baddeley
#'
#'    JmultiInhom
#'    Jdot.inhom
#'    Jcross.inhom
#'
#'    $Revision: 1.7 $ $Date: 2023/04/10 03:11:50 $

Jmulti.inhom <- function(X, I, J, 
                        lambda=NULL, lambdaI=NULL, lambdaJ=NULL,
                        lambdamin=NULL,
                        ...,
                        r=NULL, 
                        ReferenceMeasureMarkSetI=NULL,
                        ratio=FALSE){
  ## compute multitype inhomogeneous G
  ## (including determination of r and validation of lmin)
  GIJ <- GmultiInhom(X, I, J,
                     lambda, lambdaI, lambdaJ, lambdamin,
                     ...,
                     r=r,
                     ReferenceMeasureMarkSetI=ReferenceMeasureMarkSetI,
                     ratio=ratio)
  
  ## compute multitype inhomogeneous F
  FJ <- FmultiInhom(X, J,
                lambda, lambdaJ, 
                lambdamin,
                ...,
                r = GIJ$r)
  
  ## evaluate inhomogeneous J function
  if(!ratio) {
    JIJ <- eval.fv((1 - GIJ) / (1 - FJ))
  } else {
    num <- eval.fv(1 - GIJ)
    den <- eval.fv(1 - FJ)
    JIJ <- eval.fv(num / den)
    JIJ <- rat(JIJ, num, den)
  }
  
  ## relabel the fv object
  JIJ <- rebadge.fv(JIJ,
                    new.ylab  = quote(J[inhom, I, J](r)),
                    new.fname = c("J", "list(inhom,I,J)"),
                    tags      = names(JIJ),
                    new.labl  = attr(GIJ, "labl"),
                    new.yexp  = quote(J[list(inhom, I, J)](r)))
  
  ## tack on extra info
  attr(JIJ, "G") <- GIJ
  attr(JIJ, "F") <- FJ
  attr(JIJ, "dangerous") <- attr(GIJ, "dangerous")
  attr(JIJ, "conserve") <- append(attr(GIJ, "conserve"), attr(FJ, "conserve"))
  return(JIJ)
}

Jdot.inhom <- function(X, i,
                       lambdaI=NULL, lambdadot=NULL,
                       lambdamin=NULL,
                       ...,
                       r=NULL, 
                       ReferenceMeasureMarkSetI = NULL,
                       ratio = FALSE){
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE))
    stop("Point pattern must be multitype")
  marx <- marks(X, dfok = FALSE)

  if(missing(i) || is.null(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))          

  J <- rep.int(TRUE, npoints(X))
  result <- Jmulti.inhom(X, I, J, 
                         lambdaI=lambdaI, lambdaJ=lambdadot,
                         lambdamin=lambdamin,
                         ...,
                         r = r, 
                         ReferenceMeasureMarkSetI = ReferenceMeasureMarkSetI,
                         ratio = ratio)
  conserve <- attr(result, "conserve")
  result <- rebadge.as.dotfun(result, "J", "inhom", i)
  attr(result, "conserve") <- conserve
  return(result)
}

Jcross.inhom <- function(X, i, j, 
                         lambda = NULL, lambdaI = NULL, lambdaJ = NULL, 
                         lambdamin = NULL,
                         ...,
                         r = NULL, 
                         ReferenceMeasureMarkSetI = NULL,
                         ratio = FALSE) {
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE))
    stop("Point pattern must be multitype")
  marx <- marks(X, dfok=FALSE)

  if(missing(i) || is.null(i))
    i <- levels(marx)[1]
  if(missing(j) || is.null(j))
    j <- levels(marx)[2]

  I <- (marx == i)
  J <- (marx == j)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))          
  if(sum(J) == 0)
    stop(paste("No points have mark = ", j))
  
  result <- Jmulti.inhom(X, I, J, 
                         lambda, lambdaI, lambdaJ,
                         lambdamin,
                         ...,
                         r=r, 
                         ReferenceMeasureMarkSetI=ReferenceMeasureMarkSetI,
                         ratio=ratio)
  conserve <- attr(result, "conserve")
  result <- rebadge.as.crossfun(result, "J", "inhom", i, j)
  attr(result, "conserve") <- conserve
  return(result)
}

