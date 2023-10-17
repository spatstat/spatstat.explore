## clarkevans.R
## Clark-Evans statistic and test
## $Revision: 1.21 $ $Date: 2023/10/17 05:13:03 $

clarkevans <- function(X, correction=c("none", "Donnelly", "cdf"),
                       clipregion=NULL)
{
  verifyclass(X, "ppp")
  W <- X$window

  # validate correction argument
  gavecorrection <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Donnelly="Donnelly",
                             donnelly="Donnelly",
                             guard="guard",
                             cdf="cdf"),
                           multi=TRUE)

  if(("Donnelly" %in% correction) && (W$type != "rectangle")) {
    if(gavecorrection)
      warning("Donnelly correction only available for rectangular windows")
    correction <- correction[correction != "Donnelly"]
  }

  # guard correction applied iff `clipregion' is present
  isguard <- "guard" %in% correction
  askguard <- any(isguard)
  gaveguard <- !is.null(clipregion)
  if(gaveguard)
    clipregion <- as.owin(clipregion)
  if(askguard && !gaveguard) {
    warning("guard correction not performed; clipregion not specified")
    correction <- correction[!isguard]
  } else if(gaveguard && !askguard) 
    correction <- c(correction, "guard")

  result <- clarkevansCalc(X, correction, clipregion)
  if(length(result) == 1L) result <- unname(result)
  return(result)
}

clarkevans.test <- function(X, ..., 
                            correction,
                            clipregion=NULL,
                            alternative=c("two.sided", "less", "greater",
                                          "clustered", "regular"),
                            method = c("asymptotic", "MonteCarlo"),
                            nsim=999
                            ) {
  Xname <- short.deparse(substitute(X))
  miss.nsim <- missing(nsim)
  method <- match.arg(method)
  
  verifyclass(X, "ppp")
  W <- Window(X)
  nX <- npoints(X)

  if(missing(correction) || is.null(correction)) {
    correction <- switch(method,
                         MonteCarlo = "none",
                         asymptotic = if(is.rectangle(W)) "Donnelly" else "cdf")
  } else {
    #' validate SINGLE correction
    correction <- pickoption("correction", correction,
                             c(none="none",
                               Donnelly="Donnelly",
                               donnelly="Donnelly",
                               guard="guard",
                               cdf="cdf"))
  }
  
  switch(correction,
         none={
           corrblurb <- "No edge correction"
         },
         Donnelly={
           if(W$type != "rectangle")
             stop("Donnelly correction only available for rectangular windows")
           corrblurb <- "Donnelly correction"
         },
         guard={
           if(is.null(clipregion))
             stop("clipregion not specified")
           clipregion <- as.owin(clipregion)
           corrblurb <- "Guard correction"
         },
         cdf={
           corrblurb <- "CDF correction"
         })

  # alternative hypothesis
  if(missing(alternative) || is.null(alternative))
    alternative <- "two.sided"
  alternative <- pickoption("alternative", alternative,
                           c(two.sided="two.sided",
                             less="less",
                             clustered="less",
                             greater="greater",
                             regular="greater"))

  altblurb <-
    switch(alternative,
           two.sided="two-sided",
           less="clustered (R < 1)",
           greater="regular (R > 1)")

  # compute observed value
  statistic <- clarkevansCalc(X, correction=correction, clipregion=clipregion,
                              working=TRUE)
  working <- attr(statistic, "working")
  #
  switch(method,
         asymptotic = {
           #' use asymptotic standard Normal reference
           #' get appropriate standard error
           SE.R <- switch(correction,
                          none     = working[["SEnaive"]],
                          guard    = working[["SEguard"]],
                          Donnelly = working[["SEkevin"]],
                          cdf      = working[["SEcdf"]])
           #' standardised test statistic
           Z <- as.numeric((statistic - 1)/SE.R)
           p.value <- switch(alternative,
                             less=pnorm(Z),
                             greater=1 - pnorm(Z),
                             two.sided= 2*(1-pnorm(abs(Z))))
           pvblurb <- "Z-test"
         },
         MonteCarlo = {
           #' Monte Carlo p-value
           sims <- numeric(nsim)
           for(i in seq_len(nsim)) {
             Xsim <- runifpoint(nX, win=W)
             sims[i] <- clarkevansCalc(Xsim, correction=correction,
                                       clipregion=clipregion)
           }
           p.upper <- (1 + sum(sims >= statistic))/(1.0 + nsim)
           p.lower <- (1 + sum(sims <= statistic))/(1.0 + nsim)
           p.value <- switch(alternative,
                             less=p.lower,
                             greater=p.upper,
                             two.sided=min(1, 2*min(p.lower, p.upper)))
           pvblurb <- paste("Monte Carlo test based on",
                            nsim, "simulations of CSR with fixed n")
         })

  statistic <- as.numeric(statistic)
  names(statistic) <- "R"
  
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=c("Clark-Evans test", corrblurb, pvblurb),
              data.name=Xname)
  class(out) <- "htest"
  return(out)
}

clarkevansCalc <- function(X, correction="none", clipregion=NULL,
                           working=FALSE) {
  # calculations for Clark-Evans index or test
  W <- Window(X)
  areaW <- area(W)
  npts <- npoints(X)
  intensity <- npts/areaW
  # R undefined for empty point pattern
  if(npts == 0)
    return(NA)
  # Dobs = observed mean nearest neighbour distance
  nndistX <- nndist(X)
  Dobs <- mean(nndistX)
  # Dpois = Expected mean nearest neighbour distance for Poisson process
  Dpois <- 1/(2*sqrt(intensity))

  ## initialise
  statistic <- NULL
  SE.Dobs <- NULL
  if(working) {
    work <- list(areaW=areaW, npts=npts, intensity=intensity,
                 Dobs=Dobs, Dpois=Dpois)
    #' null standard error of Dobs = mean(nndist(X)) 
    SE.Dobs <- sqrt(((4-pi)*areaW)/(4 * pi))/npts  # sic
  }

  ## start computing results
  
  # Naive uncorrected value
  if("none" %in% correction) {
    Rnaive <- Dobs/Dpois
    statistic <- c(statistic, naive=Rnaive)
    if(working) {
      #' null standard error of Clark-Evans statistic Rnaive
      SE.Rnaive <- SE.Dobs / Dpois
      work <- append(work, list(SEnaive=SE.Rnaive))
    }
  }
  #' Donnelly edge correction
  if("Donnelly" %in% correction) {
     #' Edge corrected mean nearest neighbour distance, Donnelly 1978
    if(W$type == "rectangle") {
      perim <- perimeter(W)
      Dkevin  <- Dpois + (0.0514+0.0412/sqrt(npts))*perim/npts
      Rkevin <- Dobs/Dkevin
      if(working) {
        #' null standard error of adjusted Clark-Evans statistic Rkevin
        SE.Rkevin <- SE.Dobs / Dkevin
        work <- append(work,
                       list(perim=perim,
                            Dkevin=Dkevin,
                            SEkevin=SE.Rkevin))
      }
    } else {
      Rkevin <- NA
    }
    statistic <- c(statistic, Donnelly=Rkevin)
  }
  # guard area method
  if("guard" %in% correction && !is.null(clipregion)) {
    #' use nn distances from points inside `clipregion'
    ok <- inside.owin(X, , clipregion)
    Dguard <- mean(nndistX[ok])
    Rguard <- Dguard/Dpois
    statistic <- c(statistic, guard=Rguard)
    ## additional info
    if(working) {
      npts.guard <- sum(ok)
      areaWclip <- area(clipregion)
      #' null standard error of Dguard = mean(nndist(X[clipregion])) 
      SE.Dguard <- sqrt((4-pi)/(4 * pi * npts.guard * intensity))
      #' null standard error of adjusted Clark-Evans statistic Rguard
      SE.Rguard <- SE.Dguard / Dpois
      work <- append(work,
                     list(Dguard=Dguard,
                          npts.guard=npts.guard,
                          SEguard=SE.Rguard))

    }
  }
  if("cdf" %in% correction) {
    # compute mean of estimated nearest-neighbour distance distribution G
    G <- Gest(X)
    numer <- stieltjes(function(x){x}, G)$km
    denom <- stieltjes(function(x){rep.int(1, length(x))}, G)$km
    Dcdf <- numer/denom
    Rcdf <- Dcdf/Dpois
    statistic <- c(statistic, cdf=Rcdf)
    if(working) {
      #' approximate null standard error of Dobs = mean(Gest(X)) 
      SE.Dcdf <- SE.Dobs
      #' null standard error of Clark-Evans statistic Rcdf
      SE.Rcdf <- SE.Dcdf/Dpois
      work <- append(work, list(Dcdf=Dcdf,
                                SEcdf=SE.Rcdf))
    }
  }

  if(working) attr(statistic, "working") <- work

  return(statistic)
}
