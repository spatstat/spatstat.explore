\name{transect.im}
\alias{transect.im}
\title{
  Pixel Values Along a Transect
}
\description{
  Extract the pixel values of a pixel image at each point along
  a linear transect.
}
\usage{
transect.im(X, ..., from="bottomleft", to="topright",
            nsample=512, click=FALSE, add=FALSE, curve=NULL)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{from,to}{
    Optional. Start point and end point of the transect.
    Pairs of \eqn{(x,y)} coordinates
    in a format acceptable to \code{\link{xy.coords}},
    or keywords \code{"bottom"}, \code{"left"}, \code{"top"},
    \code{"right"}, \code{"bottomleft"} etc.
  }
  \item{nsample}{
    Integer. Number of sample locations along the transect.
  }
  \item{click}{
    Optional.
    Logical value.
    If \code{TRUE}, the linear transect is determined interactively
    by the user, who clicks two points on the current plot.
  }
  \item{add}{
    Logical. If \code{click=TRUE}, this argument determines
    whether to perform interactive tasks on the current plot (\code{add=TRUE})
    or to start by plotting \code{X} (\code{add=FALSE}).
  }
  \item{curve}{
    Optional. A specification of a curved transect. See the section
    on Curved Transect.
  }
}
\details{
  The pixel values of the image \code{X} along a line segment 
  will be extracted. The result is a function table (\code{"fv"} object)
  which can be plotted directly.

  If \code{click=TRUE}, then the user is prompted to click two points on
  the plot of \code{X}. These endpoints define the transect.

  Otherwise, the transect is defined by the endpoints
  \code{from} and \code{to}. The default is a diagonal transect from
  bottom left to top right of the frame.
}
\section{Curved Transect}{
  If \code{curve} is given, then the transect will be a curve.
  The argument \code{curve} should be a list with the following
  arguments:
  \describe{
    \item{f}{A function in the \R language with one argument \code{t}.}
    \item{tlim}{A numeric vector of length 2 giving the range of values
      of the argument \code{t}.}
    \item{tname}{(Optional) a character string giving the symbolic name
      of the function argument \code{t}; defaults to \code{"t"}.}
    \item{tdescrip}{(Optional) a character string giving a short
      description of the function argument \code{t}; defaults to
      \code{"curve parameter"}.}
  }
  The function \code{f} must return a 2-column matrix or data frame
  specifying the spatial coordinates \code{(x,y)} of locations along the
  curve, determined by the values of the input argument \code{t}.
}
\value{
  An object of class \code{"fv"} which can be plotted.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link[spatstat.geom]{im}}
}
\examples{
  Z <- bei.extra$elev
  plot(transect.im(Z))
}
\keyword{spatial}
\keyword{manip}
\keyword{iplot}
