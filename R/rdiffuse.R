#'
#'  rdiffuse.R
#'
#'  Random diffusion by random walk on raster
#'
#'  $Revision: 1.2 $ $Date: 2026/04/10 04:58:30 $
#'

rdiffuse <- function(X, sigma, ...) {
  UseMethod("rdiffuse")
}

rdiffuse.ppp <- function(X, sigma, ..., connect=8,
                         method=c("C", "interpreted"),
                         unround=TRUE) {
  method <- match.arg(method)
  ## >>>>. construct transition matrix etc <<<<<<<<<<<<<<<<<<<<<<
  internal <- resolve.defaults(list(setuponly=TRUE), list(...)$internal)
  stuff <- do.call(densityHeat.ppp,
                   resolve.defaults(list(x=quote(X),
                                         sigma=sigma,
                                         connect=connect,
                                         internal=internal),
                                    list(weights=NULL,
                                         show=FALSE, se=FALSE,
                                         at="pixels",
                                         leaveoneout=FALSE,
                                         extrapolate=FALSE,
                                         coarsen=FALSE),
                                    list(...)))
  ## discrete grid 
  Ximage <- stuff$Y
  ## discretised initial state of each X[i] in serialised grid
  Xpos <- stuff$Xpos
  ## transition matrix
  P <- stuff$A
  ## k-step transition matrix
  Pk <- stuff$Ak
  ## number of blocks of size k
  Nblock <- stuff$Nblock
  ## remaining number of steps
  Nrump  <- stuff$Nrump
  ## map from serial numbers 'u' to grid positions (NULL if window is rectangle)
  backmap  <- stuff$backmap
  ## convert transition matrices to efficient row-major format
  P <- as(P, "RsparseMatrix")
  Pk <- as(Pk, "RsparseMatrix")
  ## >>>>>>>>>>>>>   run parallel Markov chains <<<<<<<<<<<<<<<<<<<<<
  Z0 <- Xpos
  Zb <- runSparseMarkovChain(Pk, Z0, Nblock,
                             result="last", check=FALSE, method=method)
  Zf <- runSparseMarkovChain(P,  Zb, Nrump, 
                             result="last", check=FALSE, method=method)
  ## >>>>>>>>>>>>>  map pixel serial numbers to spatial positions <<<<<<<
  if(!is.null(backmap))
    Zf <- backmap[Zf]
  xy <- rasterxy.im(Ximage)
  xyZ <- xy[Zf,]
  ## point pattern result
  Z <- as.ppp(xyZ, W=Window(X))
  ## de-discretise
  if(unround) {
    Z <- rUnround(Z, xstep=Ximage$xstep, ystep=Ximage$ystep)
  }
  return(Z)
}

