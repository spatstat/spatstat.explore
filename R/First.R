##  spatstat.explore/R/First.R

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.explore"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatExploreVersion", vs)
  packageStartupMessage(paste("spatstat.explore", vs))
  return(invisible(NULL))
}

  
