#'
#'    rose.R
#'
#'   Rose diagrams
#'
#'   $Revision: 1.22 $  $Date: 2026/04/29 07:23:14 $
#'

rose <- function(x, ...) UseMethod("rose")

rose.default <- local({

  rose.default <- function(x, breaks = NULL, ...,
                           weights=NULL,
                           nclass=NULL,
                           unit=c("degree", "radian",
                                  "hour", "minute", "other"),
                           fullcircle=NULL, start=NULL, clockwise=NULL,
                           main) {
    if(missing(main) || is.null(main))
      main <- short.deparse(substitute(x))
    stopifnot(is.numeric(x))
    if(!is.null(weights))
      check.nvector(weights, length(x), things="observations", vname="weights")
    #' determine units and resolve auxiliary arguments
    missu <- missing(unit)
    unit <- match.arg(unit)
    unit <- validate.angles(x, unit, fullcircle=fullcircle, guess=missu)
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
                   start=start, clockwise=clockwise),
              otherargs))
  }

  graphicsAargh <- c("density", "angle", "col", "border",
                     "xlim", "ylim", "xlab", "ylab", "axes", "labels")

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
                           do.plot=TRUE, do.rings=TRUE, do.ticks=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  #' determine units
  missu <- missing(unit)
  unit <- match.arg(unit)
  #' validate and resolve arguments
  bks <- x$breaks
  unit <- validate.angles(bks, unit, fullcircle=fullcircle, guess=missu)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #' get sector sizes
  y <- x$density
  ymax <- max(y)
  #' draw disc
  insideclearance <- 0.1
  do.ticks <- !isFALSE(do.ticks) && (is.null(at) || length(at) > 0)
  outsidespace <- if(!do.ticks) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  R <- (1+insideclearance) * ymax
  DD <- disc(R)
  Rout <- (1 + outsidespace) * R
  disco <- disc(Rout)
  dont.complain.about(DD, disco)
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=quote(disco),
                                                  main=main,
                                                  type="n"), 
                                             list(...)))
  do.call.matched(plot.owin,
                  resolve.defaults(list(x=quote(DD),
                                        hatch=FALSE,
                                        add=TRUE),
                                   list(...)),
                  extrargs=graphicsPars("owin"))
  if(do.plot) {
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
      attr(result, "radii") <- rad
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
      do.call.matched(polygon, list(x=xx, y=yy, ..., col=ci))
    }
    if(do.ticks) {
      #' add tick marks
      circticks(R, at=at, unit=unit, fullcircle=fullcircle,
                start=start, clockwise=clockwise,
                labels=labels)
    }
  }
  #'
  attr(result, "histogram") <- x
  return(invisible(result))
}

rose.density <- function(x, ...,
                         unit=c("degree", "radian", "hour", "minute", "other"),
                         fullcircle=NULL,
                         start=NULL, clockwise=NULL,
                         main, labels=TRUE, at=NULL,
                         do.plot=TRUE, do.rings=TRUE, do.ticks=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- x$x
  rad <- x$y
  missu <- missing(unit)
  unit <- match.arg(unit)
  ## validate and resolve arguments
  unit <- validate.angles(ang, unit, fullcircle=fullcircle, guess=missu)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           fullcircle=fullcircle,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           do.plot=do.plot,
                           do.rings=do.rings,
                           do.ticks=do.ticks)
  return(invisible(result))
}

rose.fv <- function(x, ..., unit=c("degree", "radian",
                                   "hour", "minute", "other"),
                    fullcircle=NULL, 
                    start=NULL, clockwise=NULL,
                    main, labels=TRUE, at=NULL,
                    do.plot=TRUE, do.rings=TRUE, do.ticks=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- with(x, .x)
  rad <- with(x, .y)
  missu <- missing(unit)
  unit <- match.arg(unit)
  #' validate and resolve arguments
  unit <- validate.angles(ang, unit, fullcircle=fullcircle, guess=missu)
  a <- resolve.rose.args(unit, fullcircle=fullcircle,
                         start=start, clockwise=clockwise)
  fullcircle <- a$fullcircle
  start      <- a$start
  clockwise  <- a$clockwise
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           fullcircle=fullcircle,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           do.plot=do.plot, do.rings=do.rings,
                           do.ticks=do.ticks)
  return(invisible(result))
}

roseContinuous <- function(ang, rad, unit, ...,
                           fullcircle=NULL,
                           start=NULL, clockwise=NULL,
                           main,
                           labels=TRUE, at=NULL,
                           do.plot=TRUE, do.rings=TRUE, do.ticks=TRUE) {
  rmax <- max(rad)
  do.ticks <- !isFALSE(do.ticks) && (is.null(at) || length(at) > 0)
  #' draw disc
  insideclearance <- 0.1
  outsidespace <- if(!is.null(at) && length(at) == 0) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  R <- (1+insideclearance) * rmax
  DD <- disc(R)
  Rout <- (1 + outsidespace) * R
  disco <- disc(Rout)
  dont.complain.about(DD, disco)
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=quote(disco),
                                                  main=main,
                                                  type="n"), 
                                             list(...)))
  do.call.matched(plot.owin,
                  resolve.defaults(list(x=quote(DD),
                                        add=TRUE,
                                        hatch=FALSE),
                                   list(...)),
                  extrargs=graphicsPars("owin"),
                  skipargs="col")
  #' draw plot
  if(do.plot) {
    if(do.rings) {
      rad <- prettyinside(c(0, max(ang)), n=6)
      rad <- rad[rad > 0]
      for(radi in rad)
        plot(disc(radi), add=TRUE, border="grey", lty=2)
      attr(result, "radii") <- rad
    }
    ang <- ang2rad(ang, unit=unit, fullcircle=fullcircle,
                   start=start, clockwise=clockwise)
    xx <- rad * cos(ang)
    yy <- rad * sin(ang)
    do.call.matched(polygon, list(x=xx, y=yy, ...), extrargs="lwd")
    if(do.ticks) 
      circticks(R, at=at, unit=unit, fullcircle=fullcircle,
                start=start, clockwise=clockwise,
                labels=labels)
  }
  return(result)
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
             radian = ,
             hour = {
               at <- fullcircle * (0:23)/24
               major <- ((0:23) %% 6 == 0)
             },
             minute = {
               at <- seq(0, 55, by=15)
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

validate.angles <- function(angles,
                            unit=c("degree", "radian",
                                   "hour", "minute", "other"),
                            fullcircle=NULL, guess=TRUE) {
  #' validate
  width <- diff(range(angles))
  if(!isTRUE(guess)) {
    unit <- match.arg(unit)
  } else if(!is.null(fullcircle)) {
    unit <- "other"
  } else if(width <= 6.2832) {
    warning("Very small range of angles: treating them as radian")
    unit <- "radian"
  } else {
    unit <- "degree"
  }
  fullcircle <- resolve.rose.args(unit, fullcircle)$fullcircle
  if(width > 1.002 * fullcircle)
    stop("Range of angles exceeds a full circle")
  return(unit)
}

resolve.rose.args <- function(unit=c("degree", "radian",
                                     "hour", "minute", "other"),
                              fullcircle=NULL,
                              start=NULL, clockwise=NULL) {
  unit <- match.arg(unit)
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
  out <- list(fullcircle=fullcircle, start=start, clockwise=clockwise)
  return(out)
}
