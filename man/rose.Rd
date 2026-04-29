\name{rose}
\alias{rose}
\alias{rose.default}
\alias{rose.histogram}
\alias{rose.density}
\alias{rose.fv}
\title{Rose Diagram}
\description{
  Plots a rose diagram (rose of directions),
  the analogue of a histogram or density plot for angular data.
}
\usage{
rose(x, \dots)

\method{rose}{default}(x, breaks = NULL, \dots,
                       weights=NULL,
                       nclass = NULL,
                       unit = c("degree", "radian", "hour", "minute", "other"),
                       fullcircle = NULL, 
                       start=NULL, clockwise=NULL,
                       main)

\method{rose}{histogram}(x, \dots,
                       unit = c("degree", "radian", "hour", "minute", "other"),
                       fullcircle = NULL, 
                       start=NULL, clockwise=NULL,
                       main, col="lightgray", labels=TRUE, at=NULL,
                       do.plot = TRUE, do.rings=TRUE, do.ticks=TRUE)

\method{rose}{density}(x, \dots, 
                       unit = c("degree", "radian", "hour", "minute", "other"),
                       fullcircle = NULL, 
                       start=NULL, clockwise=NULL,
                       main, labels=TRUE, at=NULL,
                       do.plot = TRUE, do.rings=TRUE, do.ticks=TRUE)

\method{rose}{fv}(x, \dots, 
                       unit = c("degree", "radian", "hour", "minute", "other"),
                       fullcircle = NULL, 
                       start=NULL, clockwise=NULL,
                       main, labels=TRUE, at=NULL,
                       do.plot = TRUE, do.rings=TRUE, do.ticks = TRUE)
}
\arguments{
  \item{x}{
    Data to be plotted.
    A numeric vector containing angles,
    or a \code{histogram} object containing a histogram of
    angular values, or a \code{density} object containing a
    smooth density estimate for angular data,
    or an \code{fv} object giving a function
    of an angular argument.
  }
  \item{breaks, nclass}{
    Arguments passed to \code{\link[graphics]{hist}} to determine
    the histogram breakpoints. Note that \code{nclass} is the number
    of \bold{breakpoints}, so that the number of histogram bins
    is \code{nclass-1}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[graphics]{polygon}}
    controlling the appearance of the plot
    (or passed from \code{rose.default} to \code{\link[graphics]{hist}}
    to control the calculation of the histogram).
  }
  \item{unit}{
    The unit in which the angles are expressed.
    Character string (partially matched).
  }
  \item{fullcircle}{
    The period (the value of \code{x} that corresponds to one full circle)
    if \code{unit="other"}.
  }
  \item{start}{
    The starting direction for measurement of angles,
    that is, the spatial direction which corresponds to a measured angle
    of zero. Either a character string giving a compass direction
    (\code{"N"} for north, \code{"S"} for south, \code{"E"} for east,
    or \code{"W"} for west) or a number giving the angle from the
    the horizontal (East) axis to the starting direction. For example,
    if \code{unit="degree"} and \code{clockwise=FALSE},
    then \code{start=90} and \code{start="N"} are equivalent.
    The default is to measure angles anti-clockwise from the
    horizontal axis (East direction) if \code{unit="degree"} or
    \code{unit="radian"} or \code{unit="other"},
    and to measure angles clockwise from top
    (North direction) if \code{unit="hour"} or \code{unit="minute"}.
  }
  \item{clockwise}{
    Logical value indicating whether angles increase in the clockwise
    direction (\code{clockwise=TRUE}) or anti-clockwise,
    counter-clockwise direction (\code{clockwise=FALSE}).
    The default is \code{clockwise=TRUE} if
    \code{unit="hour"} or \code{unit="minute"},
    and \code{clockwise=FALSE} otherwise.
  }
  \item{weights}{
    Optional vector of numeric weights associated with \code{x}.
  }
  \item{main}{
    Optional main title for the plot.
  }
  \item{col}{
    Fill colour for the histogram polygons.
    A colour value, or a vector of colour values assigning different
    colours to each polygon (recycled to the required length).
    A colour value can be the string name of a colour
    like \code{"red"}, or an integers that refers to
    a colour in the standard palette, or a string containing
    six-letter hexadecimal codes like \code{"#F0A0FF"}.
  }
  \item{labels}{
    Either a logical value indicating whether to plot labels
    next to the tick marks, or a vector of labels for the tick marks.
  }
  \item{at}{
    Optional vector of angles at which tick marks should be plotted.
    Set \code{do.ticks=FALSE} to suppress tick marks.
  }
  \item{do.plot}{
    Logical value indicating whether to really perform the plot.
  }
  \item{do.rings}{
    Logical value indicating whether to draw concentric rings
    for scale inside the plot.
  }
  \item{do.ticks}{
    Logical value indicating whether to plot tick marks on the outer circle.
  }
}
\details{
  A rose diagram or rose of directions is the analogue of
  a histogram or bar chart for data which represent angles
  in two dimensions. The bars of the bar chart are replaced by
  circular sectors in the rose diagram.

  The function \code{rose} is generic, with a default method
  for numeric data, and methods for histograms and function tables.
  
  If \code{x} is a numeric vector, it must contain angular values
  in the range 0 to 360 (if \code{unit="degree"})
  or 0 to \code{2 * pi} (if \code{unit="radian"})
  or 0 to \code{24} (if \code{unit="hour"})
  or 0 to \code{60} (if \code{unit="minute"})
  or 0 to \code{fullcircle} (if \code{unit="other"}).
  For data expressed in hours but in a 12-hour day,
  use \code{unit="hour"} and \code{fullcircle=12}.
  A histogram of the data will first be computed (without plotting) using
  \code{\link[graphics]{hist}}. Then the rose diagram of this histogram
  will be plotted by \code{rose.histogram}.

  If \code{x} is an object of class \code{"histogram"} produced by
  the function \code{\link[graphics]{hist}}, representing the histogram
  of angular data, then the rose diagram of the densities
  (rather than the counts) in this histogram object will be plotted. 

  If \code{x} is an object of class \code{"density"} produced by
  \code{\link{circdensity}} or \code{\link[stats]{density.default}},
  representing a kernel smoothed density estimate of angular data,
  then the rose diagram of the density estimate will be plotted. 
  
  If \code{x} is a function value table (object of class \code{"fv"})
  then the argument of the function will be interpreted as an angle,
  and the value of the function will be interpreted as the radius.

  The arguments \code{start} and \code{clockwise} determine how angles
  are plotted.
  \itemize{
    \item 
    \code{start} is the spatial direction which corresponds to a measured angle
    of zero. It may be either a character string giving a compass direction
    (\code{"N"} for north, \code{"S"} for south, \code{"E"} for east,
    or \code{"W"} for west) or a number giving the angle from the
    the horizontal (East) axis to the starting direction. For example,
    if \code{unit="degree"} and \code{clockwise=FALSE},
    then \code{start=90} and \code{start="N"} are equivalent.
    \item
    \code{clockwise} is a logical value indicating
    whether angles increase in the clockwise
    direction (\code{clockwise=TRUE}) or anti-clockwise,
    counter-clockwise direction (\code{clockwise=FALSE}).
  }
  The default values of \code{start} and \code{clockwise}
  depend on \code{unit}:
  \itemize{
    \item
    If \code{unit="degree"} or \code{unit="radian"} or \code{unit="other"},
    by default the angles are interpreted using the mathematical convention
    where the zero angle is the horizontal \eqn{x} axis, and angles
    increase anti-clockwise. The default values
    are \code{start=0} and \code{clockwise=FALSE}.
    \item
    If \code{unit="hour"} or \code{unit="minute"},
    by default the angles are interpreted using the usual clock convention
    where the zero angle is the vertical \eqn{y} axis,
    and angles increase clockwise.
    The default values are \code{start="N"} and \code{clockwise=TRUE}.
  }
  Other conventions can be specified
  using the arguments \code{start} and \code{clockwise}.

  Tick marks are drawn around the outer circle if \code{do.ticks=TRUE}.
  The angular positions of the tick marks are determined
  by the argument \code{at} (in the same units as \code{x}).
  Corresponding labels are drawn at the tick marks if
  \code{labels=TRUE} (the default), and are not drawn if
  \code{labels=FALSE}. Alternatively \code{labels} can be a vector
  of labels (character strings, expressions etc) to be plotted at the
  tick marks. 
  If \code{at} is missing and \code{labels} is a vector,
  the tick marks will be equally spaced,
  with the number of tick marks equal to the length of \code{labels}.
  If \code{at} and \code{labels} are both missing,
  there is a sensible default for the tick marks and their labels.
}
\value{A window (class \code{"owin"}) containing the plotted region.}
\author{\spatstatAuthors}
\seealso{
  \code{\link{fv}}, \code{\link[graphics]{hist}},
  \code{\link{circdensity}},
  \code{\link[stats]{density.default}}.
}
\examples{
  ang <- runif(200, max=360)
  rose(ang)
  rose(ang, col=c("grey", "lightgrey"), start="N", clockwise=TRUE)

  ## hours in a 120-hour day
  hrs <- runif(200, max=12)
  Nbars <- 30
  rose(hrs, unit="h", fullcircle=12,
       labels=as.roman(c(12, 1:11)),
       nclass=Nbars+1, col=rainbow(Nbars))
}
\keyword{spatial}
\keyword{hplot}

