\name{Math.fasp}
\alias{Math.fasp}
\alias{Ops.fasp}
\alias{Complex.fasp}
\alias{Summary.fasp}
\title{S3 Group Generic Methods for Function Arrays}
\description{
  These are group generic methods for objects of class \code{"fasp"}, which
  allows for usual mathematical functions and operators to be applied
  directly to function arrays. See Details for a list of implemented functions.
}
\usage{
## S3 methods for group generics have prototypes:
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm=FALSE, drop=TRUE)}
%NAMESPACE S3method("Math", "fasp")
%NAMESPACE S3method("Ops", "fasp")
%NAMESPACE S3method("Complex", "fasp")
%NAMESPACE S3method("Summary", "fasp")
}
\arguments{
  \item{x, z, e1, e2}{objects of class \code{"fasp"}.}
  \item{\dots}{further arguments passed to methods.}
  \item{na.rm}{
    Logical value specifying whether missing values
    should be removed. 
  }
}
\details{
  Below is a list of mathematical functions and operators which are
  defined for objects of class \code{"fasp"}.
  The methods are implemented using \code{\link{eval.fasp}}, which tries
  to harmonise the functions via \code{\link{harmonise.fv}} if they aren't
  compatible to begin with.

  \enumerate{
    \item Group \code{"Math"}:
    \itemize{
      \item
      \code{abs}, \code{sign}, \code{sqrt},\cr
      \code{floor}, \code{ceiling}, \code{trunc},\cr
      \code{round}, \code{signif}

      \item
      \code{exp}, \code{log},  \code{expm1}, \code{log1p},\cr
      \code{cos}, \code{sin}, \code{tan},\cr
      \code{cospi}, \code{sinpi}, \code{tanpi},\cr
      \code{acos}, \code{asin}, \code{atan}

      \code{cosh}, \code{sinh}, \code{tanh},\cr
      \code{acosh}, \code{asinh}, \code{atanh}

      \item
      \code{lgamma}, \code{gamma}, \code{digamma}, \code{trigamma}
      \item \code{cumsum}, \code{cumprod}, \code{cummax}, \code{cummin}
    }

    \item Group \code{"Ops"}:
    \itemize{
      \item
      \code{"+"}, \code{"-"}, \code{"*"}, \code{"/"},
      \code{"^"}, \code{"\%\%"}, \code{"\%/\%"}

      \item \code{"&"}, \code{"|"}, \code{"!"}

      \item \code{"=="}, \code{"!="},
      \code{"<"}, \code{"<="}, \code{">="}, \code{">"}
    }

    \item Group \code{"Summary"}:
    \itemize{
      \item \code{all}, \code{any}
      \item \code{sum}, \code{prod}
      \item \code{min}, \code{max}
      \item \code{range}
    }

    \item Group \code{"Complex"}:
    \itemize{
      \item \code{Arg}, \code{Conj}, \code{Im}, \code{Mod}, \code{Re}
    }
  }

  For the \code{Ops} group, one of the arguments is permitted to be
  a single atomic value, or a function table, instead of a function array.
}
\seealso{
  \code{\link{eval.fasp}} for evaluating expressions involving function arrays.
}
\examples{
  ## convert array of K functions to array of L functions
  K <- alltypes(amacrine, "K")
  L <- sqrt(K/pi)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
