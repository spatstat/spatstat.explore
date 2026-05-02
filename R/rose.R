#'
#'    rose.R
#'
#'   Rose diagrams
#'
#'   $Revision: 1.32 $  $Date: 2026/05/02 05:13:29 $
#'

rose <- function(x, ...) UseMethod("rose")

rose.default <- local({

  rose.default <- function(x, breaks = NULL, ...,
                           weights=NULL,
                           nclass=NULL,
                           unit=c("degree", "radian",
                                  "hour", "minute", "other"),
                           fullcircle=NULL, start=NULL, clockwise=NULL,
                           main, add=FALSE) {
    if(missing(main) || is.null(main))
      main <- short.deparse(substitute(x))
    stopifnot(is.numeric(x))
    if(!is.null(weights))
      check.nvector(weights, length(x), things="observations", vname="weights")
    #' determine units and resolve auxiliary arguments
    missu <- missing(unit)
    unit <- resolve.angular.unit(x, unit, fullcircle=fullcircle, guess=missu)
    a <- resolve.rose.args(unit, fullcircle=fullcircle,
                           start=start, clockwise=clockwise)
    fullcircle <- a$fullcircle
    start      <- a$start
    clockwise  <- a$clockwise
    #' reduce to [0, 2pi]
    x <- x %% fullcircle
    #' determine breakpoints strictly inside full circle
    breaks <- makebreaks(x, c(0, fullcircle), breaks, nclass)
    #' histogram without weights
    h <- do.call.matched(hist.default,
                         list(x=quote(x), breaks=breaks, ..., plot=FALSE),
                         skipargs=graphicsAargh,
                         sieve=TRUE)
    result <- h$result
    otherargs <- h$otherargs
    #' redo weights, if given
    if(!is.null(weights)) {
      wh <- whist(x=x, breaks=breaks, weights=weights)
      result$count <- wh
      result$density <- wh/diff(breaks)
    }
    #
    do.call(rose.histogram,
            c(list(x=result, main=main,
                   unit=unit, fullcircle=fullcircle,
                   start=start, clockwise=clockwise, add=add),
              otherargs))
  }

  graphicsAargh <- c("density", "angle", "col", "border",
                     "xlim", "ylim", "xlab", "ylab", "axes",
                     "labels", "freq", "probability")

  makebreaks <- function(x, r, breaks=NULL, nclass=NULL) {
    use.br <- !is.null(breaks)
    if (use.br) {
      if (!is.null(nclass)) 
        warning("'nclass' is not used when 'breaks' is specified")
    } else if (!is.null(nclass) && length(nclass) == 1L) {
      breaks <- nclass
    } else breaks <- "Sturges"
    use.br <- use.br && (nB <- length(breaks)) > 1L
    if (use.br) 
      breaks <- sort(breaks)
    else {
      if (is.character(breaks)) {
        breaks <- match.arg(tolower(breaks),
                            c("sturges", 
                              "fd",
                              "freedman-diaconis",
                              "scott"))
        breaks <- switch(breaks,
                         sturges = nclass.Sturges(x), 
                         `freedman-diaconis` = ,
                         fd = nclass.FD(x),
                         scott = nclass.scott(x), 
                         stop("unknown 'breaks' algorithm"))
      }
      else if (is.function(breaks)) {
        breaks <- breaks(x)
      }
      if (length(breaks) == 1) {
        if (!is.numeric(breaks) || !is.finite(breaks) || 
            breaks < 1L) 
          stop("invalid number of 'breaks'")
        breaks <- seq(r[1], r[2], length.out=breaks)
      }
      else {
        if (!is.numeric(breaks) || length(breaks) <= 1) 
          stop(gettextf("Invalid breakpoints produced by 'breaks(x)': %s", 
                        format(breaks)), domain = NA)
        breaks <- sort(breaks)
      }
    }
    return(breaks)
  }
  
  rose.default
})


rose.histogram <- function(x, ...,
                           unit=c("degree", "radian",
                                  "hour", "minute", "other"),
                           fullcircle = NULL,
                           start=NULL, clockwise=NULL,
                           main, col="lightgray", labels=TRUE, at=NULL,
                           add=FALSE,
                           do.plot=TRUE,
                           do.circle=!add, do.rings=!add, do.ticks=!add) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ## determine type of plot: counts or probability density
  a <- resolve.hist.args(x, ...)
  freq    <- a$freq
  dotargs <- a$otherargs
  #' determine units, validate and resolve arguments
  missu <- missing(unit)
  bks <- x$breaks
  unit <- resolve.angular.unit(bks, unit, fullcircle=fullcircle, guess=missu)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #' get sector sizes
  y <- if(freq) x$counts else x$density
  ymax <- max(y)
  #' determine size of circle
  insideclearance <- 0.1
  do.ticks <- !isFALSE(do.ticks) && (is.null(at) || length(at) > 0)
  outsidespace <- if(!do.ticks) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  R <- (1+insideclearance) * ymax
  DD <- disc(R) # circle drawn outside the sectors
  Rout <- (1 + outsidespace) * R
  disco <- disc(Rout) # larger disc containing the required space
  dont.complain.about(DD, disco)
  ## create space for plot and save it as the (invisible) result
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=quote(disco),
                                                  main=main,
                                                  type="n",
                                                  do.plot=do.plot), 
                                             dotargs))
  attr(result, "histogram") <- x
  attr(result, "R") <- R
  ## 
  if(do.plot) {
    ## actually plot the circle
    if(do.circle) do.call.matched(plot.owin,
                    resolve.defaults(list(x=quote(DD),
                                          hatch=FALSE,
                                          add=TRUE),
                                     dotargs),
                    extrargs=graphicsPars("owin"))
    #' draw sectors
    ang <- ang2rad(bks, unit=unit, fullcircle=fullcircle,
                   start=start, clockwise=clockwise)
    eps <- min(abs(diff(ang)), pi/128)/2
    ## determine colours
    if(length(col) == 0) col <- "lightgray"
    nc <- length(col)
    ny <- length(y)
    if(nc != ny) {
      ## recycle colours
      col <- rep(col, ny)[seq_len(ny)]
    }
    if(do.rings) {
      rad <- prettyinside(c(0, ymax), n=6)
      rad <- rad[rad > 0]
      for(radi in rad)
        plot(disc(radi), add=TRUE, border="grey", lty=2)
      attr(result, "rings") <- rad
    }
    for(i in seq_len(ny)) {
      yi <- y[i]
      ci <- col[i]
      ## make sector coordinates
      delta <- eps * sign(ang[i+1] - ang[i])
      aa <- seq(ang[i], ang[i+1], by=delta)
      aa[length(aa)] <- ang[i+1]
      xx <- c(0, yi * cos(aa), 0)
      yy <- c(0, yi * sin(aa), 0)
      do.call.matched(polygon, c(list(x=xx, y=yy, col=ci), dotargs))
    }
    if(do.ticks) {
      #' add tick marks
      circticks(R, at=at, unit=unit, fullcircle=fullcircle,
                start=start, clockwise=clockwise,
                labels=labels)
    }
  }
  #'
  return(invisible(result))
}

rose.density <- function(x, ...,
                         unit=c("degree", "radian", "hour", "minute", "other"),
                         fullcircle=NULL,
                         start=NULL, clockwise=NULL,
                         main, labels=TRUE, at=NULL,
                         add=FALSE, do.plot=TRUE,
                         do.circle=!add, do.rings=!add, do.ticks=!add) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- x$x
  rad <- x$y
  ## validate and resolve arguments
  missu <- missing(unit)
  ci <- attr(x, "circinfo")
  unit <- resolve.angular.unit(ang,
                               unit       = unit,
                               fullcircle = fullcircle,
                               guess      = missu,
                               circinfo   = ci)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise,
                         circinfo=ci)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           fullcircle=fullcircle,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           add=add, do.plot=do.plot,
                           do.circle=do.circle,
                           do.rings=do.rings,
                           do.ticks=do.ticks)
  return(invisible(result))
}

rose.fv <- function(x, ...,
                    fmla = NULL, 
                    unit=c("degree", "radian",
                                   "hour", "minute", "other"),
                    fullcircle=NULL, 
                    start=NULL, clockwise=NULL,
                    main, labels=TRUE, at=NULL,
                    add=FALSE, do.plot=TRUE,
                    do.circle=!add, do.rings=!add, do.ticks=!add) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  #' get plot coordinates
  D <- plot(x, fmla, ..., spill=TRUE)
  ang <- D$rhsdata
  rad <- D$lhsdata # vector or matrix
  shi <- D$shind   # integer indices, length 0 or 2
  #' validate and resolve arguments
  missu <- missing(unit)
  ci <- attr(x, "circinfo")
  unit <- resolve.angular.unit(ang,
                               unit       = unit,
                               fullcircle = fullcircle,
                               guess      = missu,
                               circinfo   = ci)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise,
                         circinfo=ci)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           shind=shi,
                           fullcircle=fullcircle,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           add=add, do.plot=do.plot,
                           do.circle=do.circle, do.rings=do.rings,
                           do.ticks=do.ticks)
  return(invisible(result))
}

roseContinuous <- function(ang, rad, unit, ...,
                           fullcircle=NULL,
                           start=NULL, clockwise=NULL,
                           main,
                           labels=TRUE, at=NULL,
                           shind=NULL,
                           add=FALSE, do.plot=TRUE,
                           do.circle=!add, do.rings=!add, do.ticks=!add) {
  do.ticks <- !isFALSE(do.ticks) && (is.null(at) || length(at) > 0)
  #' determine size of circle
  insideclearance <- 0.1
  outsidespace <- if(!is.null(at) && length(at) == 0) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  rmax <- max(rad)
  R <- (1+insideclearance) * rmax
  DD <- disc(R) # circle drawn outside the curve
  Rout <- (1 + outsidespace) * R
  disco <- disc(Rout) # larger disc containing the require space
  dont.complain.about(DD, disco)
  ## create space for plot and save it as the (invisible) result
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=quote(disco),
                                                  main=main,
                                                  type="n",
                                                  do.plot=do.plot && !add), 
                                             list(...)))
  attr(result, "R") <- R
  if(do.plot) {
    ## actually plot the circle
    if(do.circle) do.call.matched(plot.owin,
                    resolve.defaults(list(x=quote(DD),
                                          add=TRUE,
                                          hatch=FALSE),
                                     list(...)),
                    extrargs=graphicsPars("owin"),
                    skipargs="col")
    if(do.rings) {
      ringrad <- prettyinside(c(0, max(rad)), n=6)
      ringrad <- ringrad[ringrad > 0]
      for(radi in ringrad)
        plot(disc(radi), add=TRUE, border="grey", lty=2)
      attr(result, "rings") <- ringrad
    }
    ang <- ang2rad(ang, unit=unit, fullcircle=fullcircle,
                   start=start, clockwise=clockwise)
    xx <- rad * cos(ang)
    yy <- rad * sin(ang)
    ncurves <- NCOL(yy)
    if(ncurves == 1) {
      do.call.matched(polygon, list(x=xx, y=yy, ...), extrargs="lwd")
    } else {
      todo <- seq_len(ncurves)
      if(length(shind) == 2) {
        ## first draw shaded region
        s1 <- shind[1]
        s2 <- shind[2]
        p1 <- list(x=xx[,s1], y=yy[,s1])
        p2 <- list(x=xx[,s2], y=yy[,s2])
        if(all(rad[,s1] <= rad[,s2])) {
          Pinside  <- p1
          Poutside <- p2
        } else {
          Pinside <- p2
          Poutside <- p1
        }
        if(Area.xypolygon(Pinside) > 0)
          Pinside <- reverse.xypolygon(Pinside)
        if(Area.xypolygon(Poutside) < 0)
          Poutside <- reverse.xypolygon(Poutside)
        SH <- owinInternalPoly(poly=list(Poutside, Pinside))
        plot(SH, add=TRUE, col="lightgrey", border="lightgrey")
        todo <- setdiff(todo, shind)
      }
      ## now draw other curves
      for(j in todo) {
        do.call.matched(lines.default,
                        list(x=xx[,j], y=yy[,j], ...),
                        extrargs=c("lwd", "col", "lty"))
      }
    }
    if(do.ticks) 
      circticks(R, at=at, unit=unit, fullcircle=fullcircle,
                start=start, clockwise=clockwise,
                labels=labels)
  }
  return(invisible(result))
}

ang2rad <- local({

  compasspoints <- c(E=0,N=90,W=180,S=270)
  
  ang2rad <- function(ang,
                      unit=c("degree", "radian", "hour", "minute", "other"),
                      fullcircle=NULL, start=NULL, clockwise=NULL) {
    ## Function to convert angular values in 'units' to radians
    ## validate arguments and resolve defaults
    unit <- match.arg(unit)
    a <- resolve.rose.args(unit, fullcircle=fullcircle,
                           start=start, clockwise=clockwise)
    fullcircle <- a$fullcircle
    start      <- a$start
    clockwise  <- a$clockwise
    ## adjust 'ang' relative to start position (expressed in 'units')
    clocksign <- if(clockwise) -1 else 1
    if(is.character(start)) {
      if(is.na(match(toupper(start), names(compasspoints))))
        stop(paste("Unrecognised compass point", sQuote(start)), call.=FALSE)
      startdegrees <- compasspoints[[start]]
      start <- switch(unit,
                      degree = startdegrees,
                      radian = pi * (startdegrees/180),
                      hour   = ,
                      minute = ,
                      other  = fullcircle * (startdegrees/360))
      # start is measured anticlockwise
      ang <- start + clocksign * ang
    } else {
      check.1.real(start)
      # start is measured according to value of 'clockwise'
      ang <- clocksign * (start + ang)
    }
    ## convert 'ang' to radians
    rad <- switch(unit,
                  degree = pi * (ang/180),
                  radian = ang,
                  hour   =,
                  minute =,
                  other  = 2 * pi * (ang/fullcircle))
    return(rad)
  }

  ang2rad
})

### ................. utilities ..................................

circticks <- function(R, at=NULL,
                      unit=c("degree", "radian", "hour", "minute", "other"),
                      fullcircle=NULL,
                      start=NULL, clockwise=NULL, labels=TRUE) {
  unit <- match.arg(unit)
  a <- resolve.rose.args(unit, fullcircle, start, clockwise)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  if(is.null(at)) {
    ## default rules for position of tick marks
    if((nlab <- length(labels)) > 1) {
      ## make 'nlab' evenly-spaced ticks
      at <- fullcircle * (0:(nlab-1))/nlab
      nat <- (at/fullcircle) * 4
      major <- abs(nat - round(nat)) < 0.01
    } else {
      ## no information - use ticks appropriate to unit
      switch(unit,
             degree = ,
             radian = {
               at <- fullcircle * (0:23)/24
               major <- ((0:23) %% 6 == 0)
             },
             hour = {
               fc <- ceiling(fullcircle)
               at <- 0:(fc-1)
               mstep <- max(primefactors(fc)) # usually 3
               major <- (at %% mstep == 0)
             },
             minute = {
               fc <- ceiling(fullcircle)
               dt <- max(primefactors(fc)) # usually 5
               at <- seq(0, fc-1, by=dt)
               major <- (at %% 15 == 0)
             },
             other = {
               at <- prettyinside(c(0, fullcircle), n=4)
               at <- at[at < fullcircle]
               nat <- (at/fullcircle) * 4
               major <- abs(nat - round(nat)) < 0.01
             })
    }
  } else {
    if(length(at) == 0) return(invisible(NULL))
    ## user specified positions of tick marks
    nat <- (at/fullcircle) * 4
    major <- abs(nat - round(nat)) < 0.01
  }
  atradians <- ang2rad(ang=at, unit=unit, fullcircle=fullcircle,
                       start=start, clockwise=clockwise)
  tx <- R * cos(atradians)
  ty <- R * sin(atradians)
  expan <- ifelse(major, 1.1, 1.05)
  segments(tx, ty, expan * tx, expan * ty, lwd=major+1)
  if(!isFALSE(labels)) {
    if(isTRUE(labels)) {
      ## default labels
      labels <- switch(unit,
                       degree = paste(round(at)),
                       radian = str2expression(
                         simplenumber(at/pi, "pi", "*", 1e-3)),
                       hour = ,
                       minute = paste(round(at)),
                       other = paste(signif(at, 3)))
    } else {
      ## user-specified labels
      if(!is.expression(labels))
        labels <- as.character(labels)
      stopifnot(length(labels) == length(at))
    }
    big <- expan + 0.1
    text(big * tx, big * ty, labels=labels)
  }
  invisible(NULL)
}

resolve.angular.unit <- function(angles,
                                 unit=c("degree", "radian",
                                        "hour", "minute", "other"),
                                 fullcircle=NULL, guess=TRUE, circinfo=NULL) {
  #' validate
  width <- diff(range(angles))
  if(!isTRUE(guess)) {
    ## unit was given 
    unit <- match.arg(unit)
  } else if(!is.null(circinfo$unit)) {
    ## from a previously calculated object
    unit <- circinfo$unit
  } else if(!is.null(fullcircle)) {
    ## unit unknown but fullcircle is given
    unit <- "other"
  } else if(width <= 6.2832) {
    ## guesswork
    warning("Very small range of angles: treating them as radian")
    unit <- "radian"
  } else {
    unit <- "degree"
  }
  return(unit)
}

resolve.rose.args <- function(unit=c("degree", "radian",
                                     "hour", "minute", "other"),
                              fullcircle=NULL,
                              start=NULL, clockwise=NULL,
                              circinfo=NULL) {
  unit <- match.arg(unit)
  fullcircle <- fullcircle %orifnull% circinfo$fullcircle
  if(unit == "other" || !is.null(fullcircle)) {
    check.1.real(fullcircle)
    stopifnot(fullcircle > 0)
  }
  if(is.null(fullcircle)) 
    fullcircle <- switch(unit,
                         degree = 360,
                         radian = 2 * pi,
                         hour   = 24,
                         minute = 60,
                         other = fullcircle)
  is.clock <- (unit %in% c("hour", "minute"))
  clockwise <- clockwise %orifnull% is.clock
  start <- start %orifnull% (if(is.clock) "N" else 0)
  out <- list(unit=unit, fullcircle=fullcircle,
              start=start, clockwise=clockwise)
  return(out)
}

resolve.hist.args <- function(x, ..., freq=NULL, probability=NULL) {
  ## imitate argument processing in hist.default
  if(!is.null(freq) && !is.null(probability) &&
     isTRUE(freq) == isTRUE(probability)) {
    warning(paste("Arguments", sQuote("freq"), "and", sQuote("probability"),
                  "were both specified, but were contradictory, so",
                  sQuote("probability"), "was ignored"),
            call.=FALSE)
  }
  freq <- !isFALSE(freq %orifnull%
                   isFALSE(probability %orifnull% isFALSE(x$equidist)))
  return(list(freq=freq, otherargs=list(...)))
}
