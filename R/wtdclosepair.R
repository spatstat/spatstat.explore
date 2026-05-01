#'
#'  wtdclosepair.R
#'
#'  $Revision: 1.2 $ $Date: 2026/05/01 02:29:02 $

weightedclosepairs <- function(X, r, correction,
                               what=c("all", "indices", "ijd")) {
  what <- match.arg(what)
  ## return list(i,j,..,weight) for all r-close pairs
  switch(correction,
         none = ,
         border = {
           cl <- closepairs(X, r, what=what)
           weight <- rep(1, length(cl$i))
         },
         isotropic = ,
         Ripley = {
           bX <- bdist.points(X)
           if(what == "indices") {
             cl <- closepairs(X, r, what="ijd")
             weight <- edge.Ripley(X[cl$i], cl$d, bdistX=bX[cl$i])
             cl <- cl[c("i", "j")]
           } else {
             cl <- closepairs(X, r, what=what)
             weight <- edge.Ripley(X[cl$i], cl$d, bdistX=bX[cl$i])
           }
         },
         translate = {
           cl <- closepairs(X, r, what="all")
           weight <- edge.Trans(dx = cl$dx,
                                dy = cl$dy,
                                W = Window(X),
                                paired=TRUE)
           switch(what,
                  indices = { cl <- cl[c("i", "j")] },
                  ijd     = { cl <- cl[c("i", "j", "d")] },
                  all     = { })
         },
         periodic = {
           cl <- closepairs(X, r, what=what, periodic=TRUE)
           weight <- rep(1, length(cl$i))
         },
         {
           warning(paste("Unrecognised correction", sQuote(correction)),
                   call.=FALSE)
           return(NULL)
         }
         )
  result <- append(cl, list(weight=as.numeric(weight)))
  return(result)
}
