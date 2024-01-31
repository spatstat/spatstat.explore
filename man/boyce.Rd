\name{boyce}
\alias{boyce}
\title{
  Boyce Index
}
\description{
  Calculate the discrete or continuous Boyce index for a spatial
  point pattern dataset.
}
\usage{
boyce(X, Z, \dots, breaks = NULL, halfwidth = NULL)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{Z}{
    Habitat suitability classes or habitat suitability index.
    Either a tessellation (object of class \code{"tess"})
    or a spatial covariate such as a pixel image (object of class
    \code{"im"}), a \code{function(x,y)} or one of the letters
    \code{"a", "b"} representing the cartesian coordinates.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{rhohat.ppp}}.
  }
  \item{breaks}{
    The breakpoint values defining discrete bands of values
    of the covariate \code{Z} for which the discrete Boyce index will
    be calculated. Either a numeric vector of breakpoints for \code{Z},
    or a single integer specifying the number of evenly-spaced
    breakpoints.
    Incompatible with \code{halfwidth}.
  }
  \item{halfwidth}{
    The half-width \eqn{h} of the interval \eqn{[z-h,z+h]} which will
    be used to calculate the continuous Boyce index \eqn{B(z)} for each
    possible value \eqn{z} of the covariate \eqn{Z}.
  }
}
\details{
  Given a spatial point pattern \code{X} and some kind of explanatory
  information \code{Z}, this function computes either the
  index originally defined by Boyce et al (2002)
  or the \sQuote{continuous Boyce index} defined by Hirzel et al (2006).

  Boyce et al (2002) defined an index of habitat suitability in which
  the study region \eqn{W} is first divided into separate subregions
  \eqn{C_1,\ldots,C_m}{C[1], C[2], ..., C[m]} based on appropriate scientific
  considerations. Then we count the number \eqn{n_j}{n[j]} of data
  points of \code{X} that fall in each subregion \eqn{C_j}{C[j]},
  measure the area \eqn{a_j}{a[j]} of each subregion \eqn{C_j}{C[j]},
  and calculate the index
  \deqn{
    B_j = \frac{n_j/n}{a_j/a}
  }{
    B[j] = (n[j]/n)/(a[j]/a)
  }
  where \eqn{a} is the total area and \eqn{n} is the total number of
  points in \code{X}.

  Hirzel et al (2006) defined another version of this index which is
  based on a continuous spatial covariate. For each possible value \eqn{z}
  of the covariate \eqn{Z}, 
  consider the region \eqn{C(z)} where the value of the covariate
  lies between \eqn{z-h} and \eqn{z+h}, where \eqn{h} is the
  chosen \sQuote{halfwidth}. The \sQuote{continuous Boyce index} is
  \deqn{
    B(z) = \frac{n(z)/n}{a(z)/a}
  }{
    B(z) = (n(z)/n)/(a(z)/a)
  }
  where \eqn{n(z)} is the number of points of \code{X}
  falling in \eqn{C(z)}, and \eqn{a(z)} is  the area of \eqn{C(z)}.

  If \code{Z} is a tessellation (object of class \code{"tess"}),
  the algorithm calculates the original (\sQuote{discrete}) Boyce index
  (Boyce et al, 2002)
  for each tile of the tessellation. The result is another tessellation,
  identical to \code{Z} except that the mark values are the
  values of the discrete Boyce index.

  If \code{Z} is a pixel image whose values are categorical (i.e. factor
  values), then \code{Z} is treated as a tessellation, with one tile
  for each level of the factor. The discrete Boyce index is then
  calculated. The result is a tessellation with marks that are the
  values of the discrete Boyce index.
  
  Otherwise, if \code{Z} is a spatial covariate such as a pixel image,
  a \code{function(x,y)} or one of the characters \code{"x"} or
  \code{"y"}, then exactly one of the arguments \code{breaks} or
  \code{halfwidth} must be given. 
  \itemize{
    \item if \code{halfwidth} is given, it should be a single positive
    number. The continuous Boyce index (Hirzel et al, 2006)
    is computed using the specified halfwidth \eqn{h}.
    The result is an object of class \code{"fv"} that can be plotted
    to show \eqn{B(z)} as a function of \eqn{z}.
    \item if \code{breaks} is given, it can be either a numeric vector
    of possible values of \code{Z} defining the breakpoints for the
    bands of values of \code{Z}, or a single integer specifying the
    number of evenly-spaced breakpoints that should be created.
    The discrete Boyce index (Boyce et al, 2002) is computed.
    The result is an object of class \code{"fv"} that can be plotted
    to show the discrete Boyce index as a function of \eqn{z}.
  }
  When \code{Z} is a spatial covariate (not factor-valued), the calculation is performed
  using \code{\link{rhohat.ppp}} (since the Boyce index is a special case
  of \code{rhohat}). Arguments \code{\dots} passed to
  \code{\link{rhohat.ppp}} control the accuracy of the spatial discretisation
  and other parameters of the algorithm.
}
\value{
  A tessellation (object of class \code{"tess"})
  or a function value table (object of class \code{"fv"})
  as explained above.
}
\references{
Boyce, M.S., Vernier, P.R., Nielsen, S.E. and Schmiegelow, F.K.A. (2002)
Evaluating resource selection functions.
\emph{Ecological modelling} \bold{157}, 281--300.

Hirzel, A.H., Le Lay, V., Helfer, V., Randin, C. and Guisan, A. (2006)
Evaluating the ability of habitat suitability models
to predict species presences.
\emph{Ecological Modelling} \bold{199}, 142--152.
}
\author{
  \adrian
}
\seealso{
  \code{\link{rhohat}}
}
\examples{
  online <- interactive()
  ## a simple tessellation
  V <- quadrats(Window(bei), 4, 3)
  if(online) plot(V)

  ## discrete Boyce index for a simple tessellation
  A <- boyce(bei, V)

  if(online) {
   plot(A, do.col=TRUE)
   marks(A)
   tilenames(A)
  }

  ## spatial covariate: terrain elevation
  Z <- bei.extra$elev

  ## continuous Boyce index for terrain elevation
  BC <- boyce(bei, Z, halfwidth=10)

  if(online) plot(BC)

  ## discrete Boyce index for terrain elevation steps of height 5 metres
  bk <- c(seq(min(Z), max(Z), by=5), Inf)
  BD <- boyce(bei, Z, breaks=bk)

  if(online) plot(BD)
}
\keyword{spatial}
\keyword{nonparametric}
