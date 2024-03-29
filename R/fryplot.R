#
#  fryplot.R
#
#  $Revision: 1.19 $ $Date: 2024/02/04 08:04:51 $
#

fryplot <- function(X, ..., width=NULL, from=NULL, to=NULL, axes=FALSE) {
  Xname <- short.deparse(substitute(X))
  X <- as.ppp(X)
  b <- as.rectangle(X)
  halfspan <- with(b, c(diff(xrange), diff(yrange)))
  if(!is.null(width)) {
    halfwidth <- ensure2vector(width)/2
    halfspan <- pmin.int(halfspan, halfwidth)
  }
  bb <- owinInternalRect(c(-1,1) * halfspan[1L], c(-1,1) * halfspan[2L])
  Y <- frypoints(X, from=from, to=to, dmax=diameter(bb))[bb]
  dont.complain.about(Y)
  do.call(plot.ppp,
          resolve.defaults(list(x=quote(Y)),
                           list(...),
                           list(main=paste("Fry plot of", Xname))))
  if(axes) {
    lines(c(0,0), c(-1,1) * halfspan[1L])
    lines(c(-1,1) * halfspan[2L], c(0,0))
  }
  return(invisible(NULL))
}

frypoints <- function(X, from=NULL, to=NULL, dmax=Inf) {
  X <- as.ppp(X)
  b <- as.rectangle(X)
  bb <- owinInternalRect(c(-1,1) * diff(b$xrange), c(-1,1) * diff(b$yrange))
  n <- X$n
  xx <- X$x
  yy <- X$y
  ## determine (dx, dy) for all relevant pairs
  if(is.null(from) && is.null(to)) {
    if(is.infinite(dmax)) {
      dx <- outer(xx, xx, "-")
      dy <- outer(yy, yy, "-")
      notsame <- matrix(TRUE, n, n)
      diag(notsame) <- FALSE
      DX <- as.vector(dx[notsame])
      DY <- as.vector(dy[notsame])
      I <- row(notsame)[notsame]
      J <- col(notsame)[notsame]
    } else {
      cl <- closepairs(X, dmax)
      DX <- cl$dx
      DY <- cl$dy
      I  <- cl$j  ## sic: I is the index of the 'TO' element
      J  <- cl$i  ## sic
    }
  } else {
    seqn <- seq_len(n)
    from <- if(is.null(from)) seqn else seqn[from]
    to   <- if(is.null(to))   seqn else seqn[to]
    if(is.infinite(dmax)) {
      dx <- outer(xx[to], xx[from], "-")
      dy <- outer(yy[to], yy[from], "-")
      notsame <- matrix(TRUE, n, n)
      diag(notsame) <- FALSE
      notsame <- notsame[to, from, drop=FALSE]
      DX <- as.vector(dx[notsame])
      DY <- as.vector(dy[notsame])
      I <- row(notsame)[notsame]
      J <- col(notsame)[notsame]
    } else {
      cl <- crosspairs(X[from], X[to], dmax)
      ok <- with(cl, from[i] != to[j])
      DX <- cl$dx[ok]
      DY <- cl$dy[ok]
      I  <- cl$j[ok]
      J  <- cl$i[ok]
    }
  }
  ## form into point pattern
  Fry <- ppp(DX, DY, window=bb, check=FALSE)
  if(is.marked(X)) {
    marx <- as.data.frame(marks(X))
    marxto <- if(is.null(to)) marx else marx[to, ,drop=FALSE]
    marks(Fry) <- marxto[I, ]
  }
  attr(Fry, "indices") <- list(I=I, J=J)
  return(Fry)
}
