#'
#'   bw.CvL.adaptive.R
#'
#'   $Revision: 1.8 $ $Date: 2022/06/25 04:31:57 $
#'
#'   Original code by Marie-Colette van Lieshout
#'   Modified by Adrian Baddeley
#' 
#'   Copyright (c) Marie-Colette van Lieshout and Adrian Baddeley 2022
#'   GNU Public Licence >= 2.0

bw.CvL.adaptive <- function(X, ..., 
                            hrange=NULL, nh=16, h=NULL,
                            bwPilot=bw.scott.iso(X),
                            edge=FALSE, diggle=TRUE)
{
   verifyclass(X, "ppp")

   W <- Window(X)
   lW <- area.owin(W)

   if(!is.null(h)) {
     stopifnot(is.numeric(h))
     stopifnot(all(h > 0))
   } else {
     ## determine range of h
     if(!is.null(hrange)) {
       check.range(hrange)
       if(any(hrange <= 0)) stop("All h values must be positive")
     } else {
       nnd <- nndist(X)
       hrange <- c(min(nnd[nnd > 0]), diameter(as.owin(X))/2)
     }
     check.1.integer(nh)
     stopifnot(nh > 1)
     h <- geomseq(from=hrange[1L], to=hrange[2L], length.out=nh)
   }

   if(!is.null(bwPilot)) {
     check.1.real(bwPilot)
     stopifnot(bwPilot > 0)
   }
      
   pdens <- density(X, sigma=bwPilot,
                    edge=TRUE, diggle=TRUE, at="pixels", leaveoneout=FALSE)

   lp2 <- cv <- numeric(nh)
   for (i in 1:nh) {
     lamxi <- adaptive.density(X, h0=h[i],
                               pilot=pdens,
                               method="kernel", 
                               edge=edge, diggle=diggle, at = "points",
                               leaveoneout = FALSE,
                               ...)
     cv[i] <- sum(1/lamxi)
     lp2[i] <- (cv[i] - lW)^2
   }
   result <- bw.optim(lp2, h,
                      optimum="min",
                      cvname="lp2", hname="h",
                      criterion="Cronie-Van Lieshout",
                      unitname=unitname(X),
                      CvL=cv)
   return(result)
}
