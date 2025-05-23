\name{plot.fv}
\alias{plot.fv}
\title{Plot Function Values}
\description{
  Plot method for the class \code{"fv"}.
}
\usage{
 \method{plot}{fv}(x, fmla, \dots, subset=NULL, lty=NULL, col=NULL, lwd=NULL,
           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
           clip.xlim=TRUE, ylim.covers=NULL,
           legend=!add, legendpos="topleft", legendavoid=missing(legendpos),
           legendmath=TRUE, legendargs=list(),
           shade=fvnames(x, ".s"), shadecol="grey",
           add=FALSE, log="",
           mathfont=c("italic", "plain", "bold", "bolditalic"), 
           limitsonly=FALSE, do.plot=TRUE)
}
\arguments{
  \item{x}{
    An object of class \code{"fv"}, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{fmla}{
    an R language formula 
    determining which variables or expressions are plotted.
    Either a formula object, or a string that can be parsed as a
    formula.
    See Details.
  }
  \item{subset}{
    (optional) subset of rows of the data frame that will be plotted.
  }
  \item{lty}{
    (optional) numeric vector of values of the graphical parameter
    \code{lty} controlling the line style of each plot.
  }
  \item{col}{
    (optional) numeric vector of values of the graphical parameter
    \code{col} controlling the colour of each plot.
  }
  \item{lwd}{
    (optional) numeric vector of values of the graphical parameter
    \code{lwd} controlling the line width of each plot.
  }
  \item{xlim}{
    (optional) range of x axis
  }
  \item{ylim}{
    (optional) range of y axis
  }
  \item{xlab}{
    (optional) label for x axis
  }
  \item{ylab}{
    (optional) label for y axis
  }
  \item{\dots}{
    Extra arguments passed to \code{plot.default}.
  }
  \item{clip.xlim}{
    Logical value specifying whether the range of the horizontal axis
    \code{xlim} should be automatically restricted to a subset of the
    range of the available data. See the section on
    \bold{Controlling the horizontal axis limits} below.
  }
  \item{ylim.covers}{
    Optional vector of \eqn{y} values that must be included in the
    \eqn{y} axis. For example \code{ylim.covers=0} will ensure that the
    \eqn{y} axis includes the origin.
  }
  \item{legend}{
    Logical flag or \code{NULL}. If \code{legend=TRUE}, the algorithm
    plots a legend in the top left corner of the plot,
    explaining the meaning of the different line types and colours.
  }
  \item{legendpos}{
    The position of the legend. Either a character string keyword
    (see \code{\link[graphics]{legend}} for keyword options)
    or a pair of coordinates in the format \code{list(x,y)}.
    Alternatively if \code{legendpos="float"}, a location will be
    selected inside the plot region, avoiding the graphics.
  }
  \item{legendavoid}{
    Whether to avoid collisions between the legend and the graphics.
    Logical value.
    If \code{TRUE}, the code will check for collisions between the
    legend box and the graphics, and will override \code{legendpos}
    if a collision occurs.
    If \code{FALSE}, the value of \code{legendpos} is always respected.
  }
  \item{legendmath}{
    Logical. If \code{TRUE}, the legend will display the
    mathematical notation for each curve. If \code{FALSE}, the legend text
    is the identifier (column name) for each curve.
  }
  \item{legendargs}{
    Named list containing additional arguments to be passed to
    \code{\link{legend}} controlling the appearance of the legend.
  }
  \item{shade}{
    A character vector giving the names of two columns of \code{x},
    or another type of index that identifies two columns.
    When the corresponding curves are plotted,
    the region between the curves will be shaded in light grey.
    The object \code{x} may or may not contain two columns which are designated
    as boundaries for shading; they are identified by \code{fvnames(x, ".s")}.
    The default is to shade between these two curves
    if they exist. To suppress this behaviour, set \code{shade=NULL}.
  }
  \item{shadecol}{
    The colour to be used in the \code{shade} plot.
    A character string or an integer specifying a colour.
  }
  \item{add}{
    Logical. Whether the plot should be added to an existing plot
  }
  \item{log}{
    A character string which contains \code{"x"} if the x axis is to
    be logarithmic, \code{"y"} if the y axis is to be logarithmic and
    \code{"xy"} or \code{"yx"} if both axes are to be logarithmic.
  }
  \item{mathfont}{
    Character string. The font to be used for mathematical
    expressions in the axis labels and the legend.
  }
  \item{limitsonly}{
    Logical. If \code{FALSE}, plotting is performed normally.
    If \code{TRUE}, no plotting is performed at all; 
    just the \eqn{x} and \eqn{y} limits of the plot are computed
    and returned. 
  }
  \item{do.plot}{
    Logical value indicating whether to actually plot the graph of \code{x}.
    Setting \code{do.plot=FALSE} will simply return the data frame
    giving the meaning of the different line types and colours,
    without plotting them.
  }
}
\value{
  Invisible: either \code{NULL}, or a data frame giving the meaning of the
  different line types and colours.
}
\details{
  This is the \code{plot} method for the class \code{"fv"}.

  An object of class \code{"fv"} is a convenient way of storing several different
  statistical estimates of a summary function; see \code{\link{fv.object}}.
  The default behaviour, executed by \code{plot(x)}, displays these
  different estimates as curves with different colours and line styles,
  and plots a legend explaining them.

  The use of the argument \code{fmla} is like \code{plot.formula}, but offers
  some extra functionality.

  The left and right hand sides of \code{fmla} are evaluated,
  and the results are plotted against each other
  (the left side on the \eqn{y} axis 
  against the right side on the \eqn{x} axis).

  The left and right hand sides of \code{fmla} may be
  the names of columns of the data frame \code{x},
  or expressions involving these names. If a variable in \code{fmla}
  is not the name of a column of \code{x}, the algorithm will search for
  an object of this name in the environment where \code{plot.fv} was
  called, and then in the enclosing environment, and so on. 

  Multiple curves may be specified by a single formula
  of the form 
  \code{cbind(y1,y2,\dots,yn) ~ x}, where \code{x,y1,y2,\dots,yn} are
  expressions involving the variables in the data frame.
  Each of the variables \code{y1,y2,\dots,yn} in turn will be plotted
  against \code{x}. 
  See the examples.

  Convenient abbreviations which can be used in the formula
  are 
  \itemize{
    \item the symbol \code{.} which represents all the
    columns in the data frame that will be plotted by default;
    \item the symbol \code{.x} which represents the function argument;
    \item the symbol \code{.y} which represents the recommended value
    of the function.
  }
  For further information, see \code{\link{fvnames}}.

  The value returned by this plot function indicates the
  meaning of the line types and colours in the plot. It can be used
  to make a suitable legend for the plot if you want to do this
  by hand. See the examples.

  The argument \code{shade} can be used to display critical bands
  or confidence intervals. If it is not \code{NULL}, then it should be
  a subset index for the columns of \code{x}, that identifies exactly
  2 columns. When the corresponding curves are plotted, the region
  between the curves will be shaded in light grey. See the Examples.

  The default values of \code{lty}, \code{col} and \code{lwd} can
  be changed using \code{\link[spatstat.geom]{spatstat.options}("plot.fv")}.

  Use \code{type = "n"} to create the plot region and draw the axes
  without plotting any curves.

  Use \code{do.plot=FALSE} to suppress all plotting.
  The return value is a data frame giving the meaning of the
  different line types and colours which would have been plotted.

  Use \code{limitsonly=TRUE} to suppress all plotting
  and just compute the \eqn{x} and \eqn{y} limits. This can be used
  to calculate common \eqn{x} and \eqn{y} scales for several plots.

  To change the kind of parenthesis enclosing the
  explanatory text about the unit of length, use
  \code{\link[spatstat.geom]{spatstat.options}('units.paren')}
}
\section{Controlling the horizontal axis limits}{
  The plot generated by \code{plot(x)}
  does not necessarily display all the data that is contained in the object.
  The range of values of the function argument \eqn{r}
  displayed in the plot may be narrower than the 
  range of values actually contained in the data frame.

  To override this behaviour and display all the available data,
  set \code{clip.xlim=FALSE}.

  Statistical literature for summary functions of spatial data
  recommends that, when the function is plotted,
  the values of the function argument on the horizontal axis
  should be restricted to a limited range of values.
  For example, Ripley recommends that the K-function \eqn{K(r)}
  should be plotted only for
  values of distance \eqn{r} between \eqn{0} and \eqn{b/4} where \eqn{b}
  is the shortest side of the enclosing rectangle of the data.

  This may be desirable so that the interesting detail is clearly visible in the
  plot. 
  It may be necessary because values outside the recommended range
  are theoretically invalid, or unreliable due to high variance or large bias.
  
  To support this standard practice, each object of class \code{"fv"}
  may include data specifying a \dQuote{recommended range} of values of 
  the function argument. The object produced by \code{\link{Kest}}
  includes a recommended range following Ripley's recommendation above.
  Similarly for \code{\link{Gest}}, \code{\link{Fest}} and many other
  commands.
  
  When \code{plot(x)} is executed, the horizontal axis is restricted
  to the recommended range of values. This recommendation can be
  overridden by setting \code{clip.xlim=FALSE}
  or by specifying the numerical limits \code{xlim}.
}
\examples{
   K <- Kest(cells)
   # K is an object of class "fv"

   plot(K, iso ~ r)                # plots iso against r

   plot(K, sqrt(iso/pi) ~ r)   # plots sqrt(iso/r)  against r

   plot(K, cbind(iso,theo) ~ r)   # plots iso against r  AND theo against r

   plot(K, .  ~ r)            # plots all available estimates of K against r

   plot(K, sqrt(./pi) ~ r)   # plots all estimates of L-function
                             # L(r) = sqrt(K(r)/pi)

   plot(K, cbind(iso,theo) ~ r, col=c(2,3))
                                   # plots iso against r  in colour 2
                                   # and theo against r in colour 3

   plot(K, iso ~ r, subset=quote(r < 0.2))
                                   # plots iso against r for r < 10

   # Can't remember the names of the columns? No problem..
   plot(K, sqrt(./pi) ~ .x)

   # making a legend by hand
   v <- plot(K, . ~ r, legend=FALSE)
   legend("topleft", legend=v$meaning, lty=v$lty, col=v$col)

   # significance bands
   KE <- envelope(cells, Kest, nsim=19)
   plot(KE, shade=c("hi", "lo"))

   # how to display two functions on a common scale
   Kr <- Kest(redwood)
   a <- plot(K, limitsonly=TRUE)
   b <- plot(Kr, limitsonly=TRUE)
   xlim <- range(a$xlim, b$xlim)
   ylim <- range(a$ylim, b$ylim)
   opa <- par(mfrow=c(1,2))
   plot(K, xlim=xlim, ylim=ylim)
   plot(Kr, xlim=xlim, ylim=ylim)
   par(opa)
   # For a shortcut, try plot(anylist(K, Kr), equal.scales=TRUE)
}
\seealso{
  \code{\link{fv.object}},
  \code{\link[spatstat.explore]{Kest}}
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{hplot}
 
