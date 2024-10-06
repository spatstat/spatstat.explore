#'
#' blurHeat.R
#' 
#'   Image blurring by diffusion
#' 
#'   Copyright (C) 2018-2024 Adrian Baddeley, Tilman Davies and Suman Rakshit
#' 
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.3 $ $Date: 2024/10/06 02:28:55 $

blurHeat <- function(X, ...) {
  UseMethod("blurHeat")
}

blurHeat.im <- function(X, sigma, ..., connect=8,
                        symmetric=FALSE, k=1, show=FALSE) {
  Y <- as.im(X)
  check.1.integer(k)
  stopifnot(k >= 1)
  if(!(connect %in% c(4,8)))
    stop("connectivity must be 4 or 8")
  if(is.im(sigma)) {
    # ensure Y and sigma are on the same grid
    A <- harmonise(Y=Y, sigma=sigma)
    Y <- A$Y
    sigma <- A$sigma
  } else if(is.function(sigma)) {
    sigma <- as.im(sigma, as.owin(Y))
  } else check.1.real(sigma)
  #' initial state
  v <- as.matrix(Y)
  u <- as.vector(v)
  #' symmetric random walk?
  if(symmetric) {
    asprat <- with(Y, ystep/xstep)
    if(abs(asprat-1) > 0.01)
      warning(paste("Symmetric random walk on a non-square grid",
                    paren(paste("aspect ratio", asprat))),
              call.=FALSE)
  }
  #' determine appropriate jump probabilities & time step
  pmax <- 1/(connect+1) # maximum permitted jump probability
  xstep <- Y$xstep
  ystep <- Y$ystep
  minstep <- min(xstep, ystep)
  if(symmetric) {
    #' all permissible transitions have the same probability 'pjump'.
    #' Determine Nstep, and dt=sigma^2/Nstep, such that
    #' Nstep >= 16 and M * pjump * minstep^2 = dt
    M <- if(connect == 4) 2 else 6
    Nstep <- max(16, ceiling(max(sigma)^2/(M * pmax * minstep^2)))    
    sn <- (sigma^2)/Nstep
    px <- py <- pxy <- sn/(M * minstep^2)
  } else {
    #' px is the probability of jumping 1 step to the right
    #' py is the probability of jumping 1 step up
    #' if connect=4, horizontal and vertical jumps are exclusive.
    #' if connect=8, horizontal and vertical increments are independent
    #' Determine Nstep, and dt = sigma^2/Nstep, such that
    #' Nstep >= 16 and 2 * pmax * minstep^2 = dt
    Nstep <- max(16, ceiling(max(sigma)^2/(2 * pmax * minstep^2)))
    sn <- (sigma^2)/Nstep
    px <- sn/(2 * xstep^2)
    py <- sn/(2 * ystep^2)
    if(max(px) > pmax) stop("Internal error: px exceeds pmax")
    if(max(py) > pmax) stop("Internal error: py exceeds pmax")
    if(connect == 8) pxy <- px * py
  }
  #' construct adjacency matrices
  dimv <- dim(v)
  my <- gridadjacencymatrix(dimv, across=FALSE, down=TRUE, diagonal=FALSE)
  mx <- gridadjacencymatrix(dimv, across=TRUE,  down=FALSE, diagonal=FALSE)
  if(connect == 8)
    mxy <- gridadjacencymatrix(dimv, across=FALSE,  down=FALSE, diagonal=TRUE)
  #' restrict to window
  if(anyNA(u)) {
    ok <- !is.na(u)
    u <- u[ok]
    mx <- mx[ok,ok,drop=FALSE]
    my <- my[ok,ok,drop=FALSE]
    if(connect == 8) 
      mxy <- mxy[ok,ok,drop=FALSE]
    if(is.im(sigma)) {
      px <- px[ok]
      py <- py[ok]
      if(connect == 8) 
        pxy <- pxy[ok]
    }
  } else ok <- TRUE
  #' construct iteration matrix
  if(connect == 4) {
    A <- px * mx + py * my
  } else {
    A <- px * (1 - 2 * py) * mx + py * (1 - 2 * px) * my + pxy * mxy
  }
  #' debug
  stopifnot(min(rowSums(A)) >= 0)
  stopifnot(max(rowSums(A)) <= 1)
  #' 
  diag(A) <- 1 - rowSums(A)
  #' k-step transition probabilities
  if(k > 1) {
    Ak <- A
    for(j in 2:k) Ak <- Ak %*% A
  } else Ak <- A
  k <- as.integer(k)
  Nstep <- as.integer(Nstep)
  Nblock <- Nstep/k
  Nrump  <- Nstep - Nblock * k
  #' run
  U <- u
  Z <- Y
  if(!show) {
    for(istep in 1:Nblock) U <- U %*% Ak
  } else {
    opa <- par(ask=FALSE)
    each <- max(1, round(Nblock/60))
    for(istep in 1:Nblock) {
      U <- U %*% Ak
      if(istep %% each == 0) {
        Z[] <- as.vector(U)
        f <- sqrt(istep/Nstep)
        main <- if(is.im(sigma)) paste(signif(f, 3), "* sigma") else
                paste("sigma =", signif(f * sigma, 3))
        plot(Z, main=main)
        Sys.sleep(0.4)
      }
    }
    par(opa)
  }
  if(Nrump > 0) for(istep in 1:Nrump) U <- U %*% A
  #' pack up
  Z[] <- as.vector(U)
  return(Z)
}

