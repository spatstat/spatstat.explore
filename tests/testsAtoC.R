#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.explore
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.explore)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
## tests/cdf.test.R


local({
  NSIM <- 9
  op <- spatstat.options(ndummy.min=16, npixel=32)
  AA <- split(ants, un=FALSE)
  AC <- AA[["Cataglyphis"]]
  AM <- AA[["Messor"]]
  DM <- distmap(AM)
  if(ALWAYS) {
    ## Check cdf.test with strange data
    ## Marked point patterns with some marks not represented
    ## should produce a warning, rather than a crash:
    cdf.test(AC, DM)
  }
  if(FULLTEST) {
    ## should be OK:
    cdf.test(unmark(AC), DM)
    cdf.test(unmark(AC), DM, "cvm")
    cdf.test(unmark(AC), DM, "ad")
    ## other code blocks
    cdf.test(finpines, "x")
  }
})


#'    tests/circular.R
#'
#'    Circular data and periodic distributions
#'
#'    $Revision: 1.4 $  $Date: 2020/04/28 12:58:26 $


local({
  if(ALWAYS) {
    a <- pairorient(redwood, 0.05, 0.15, correction="none")
    rose(a)
  }
  if(FULLTEST) {
    b <- pairorient(redwood, 0.05, 0.15, correction="best")
    rose(b, start="N", clockwise=TRUE)
  }
  if(ALWAYS) {
    #' arcs on the circle 
    #'       (depends on numerical behaviour)
    set.seed(19171025)
    aa <- replicate(7, runif(1, 0, 2*pi) + c(0, runif(1, 0, pi)),
                    simplify=FALSE)
    bb <- circunion(aa)

    assertsingle <- function(x, a, id) {
      y <- circunion(x)
      if(length(y) != 1 || max(abs(y[[1]] - a)) > .Machine$double.eps)
        stop(paste("Incorrect result from circunion in case", id),
             call.=FALSE)
      invisible(NULL)
    }

    assertsingle(list(c(pi/3, pi), c(pi/2, 3*pi/2)),
                 c(pi/3, 3*pi/2),
                 1)
    assertsingle(list(c(0, pi/2), c(pi/4, pi)),
                 c(0,pi),
                 2)
    assertsingle(list(c(-pi/4, pi/2), c(pi/4, pi)),
                 c((2-1/4)*pi, pi),
                 3)
  }
})

  
#'
#'   tests/closecore.R
#'
#' check 'closepairs/crosspairs' code
#' invoked in core package
#'
#' $Revision: 1.4 $ $Date: 2021/04/17 04:16:43 $
#' 
#' ------- All this code must be run on every hardware -------
#'

local({
  #' weightedclosepairs is in wtdclosepair.R
  wi <- weightedclosepairs(redwood, 0.05, "isotropic")
  if(FULLTEST) {
    wt <- weightedclosepairs(redwood, 0.05, "translate")
    wp <- weightedclosepairs(redwood, 0.05, "periodic")
  }
  #' markmarkscatter uses closepairs.pp3
  X <- runifpoint3(100)
  marks(X) <- runif(100)
  markmarkscatter(X, 0.2)
  if(FULLTEST) {
    markmarkscatter(X[FALSE], 0.2)
  }
})

#'
#'     contact.R
#'
#'   Check machinery for first contact distributions
#'
#'   $Revision: 1.8 $  $Date: 2021/04/17 02:25:55 $

local({
  if(ALWAYS) {
    #' reduce complexity
    Y <- as.mask(heather$coarse, dimyx=c(50, 25))
    
    X <- runifpoint(100, win = complement.owin(Y))
    if(FULLTEST) G <- Gfox(X, Y)
    J <- Jfox(X, Y)

    Y <- as.polygonal(Y)
    X <- runifpoint(100, win = complement.owin(Y))
    if(FULLTEST) G <- Gfox(X, Y)
    J <- Jfox(X, Y)

    op <- spatstat.options(exactdt.checks.data=TRUE)
    U <- exactdt(X)
    spatstat.options(op)
  }
})

reset.spatstat.options()
