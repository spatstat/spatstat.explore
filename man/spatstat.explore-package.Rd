\name{spatstat.explore-package}
\alias{spatstat.explore-package} 
\alias{spatstat.explore} 
\docType{package}
\title{The spatstat.explore Package}
\description{
  The \pkg{spatstat.explore} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the core functionality
  for statistical analysis and modelling of spatial data.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.

  This sub-package \pkg{spatstat.explore} contains the 
  user-level functions
  that perform exploratory data analysis and
  nonparametric data analysis of spatial data.

  (The main exception is that functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.)
}
\section{Structure of the spatstat family}{  
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.univar} containing functions for estimating
    probability distributions of random variables
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.explore} containing the functionality
    for exploratory data analysis and nonparametric analysis of spatial data.
    \item \pkg{spatstat.model} containing the functionality
    for statistical modelling, model-fitting, formal statistical
    inference  and informal model diagnostics.
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.

  For an overview of all the functions available in
  the sub-packages of \pkg{spatstat},
  see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package.

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Overview of Functionality in \pkg{spatstat.explore}}{
  The \pkg{spatstat} family of packages
  is designed to support a complete statistical analysis
  of spatial data. It supports
  \itemize{
    \item creation, manipulation and plotting of point patterns;
    \item exploratory data analysis;
    \item spatial random sampling;
    \item simulation of point process models;
    \item parametric model-fitting;
    \item non-parametric smoothing and regression;
    \item formal inference (hypothesis tests, confidence intervals);
    \item model diagnostics.
  }
  For an overview, see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package. 

  Following is a list of the functionality provided in the
  \pkg{spatstat.explore} package only.

  \bold{To simulate a random point pattern:}
  
   Functions for generating random point patterns
   are now contained in the \pkg{spatstat.random} package.

      \bold{To interrogate a point pattern:}

      \tabular{ll}{
	\code{\link[spatstat.explore]{density.ppp}} \tab
	kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.explore]{densityHeat.ppp}} \tab
	diffusion kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.explore]{Smooth.ppp}} \tab
	kernel smoothing of marks of point pattern\cr
	\code{\link[spatstat.explore]{sharpen.ppp}} \tab
	data sharpening\cr
      }

      \bold{Manipulation of pixel images:}

      An object of class \code{"im"} represents a pixel image.

      \tabular{ll}{
	\code{\link[spatstat.explore]{blur}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.explore]{Smooth.im}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.explore]{transect.im}} \tab line transect of image \cr
	\code{\link[spatstat.geom]{pixelcentres}} \tab extract centres of pixels \cr
	\code{\link[spatstat.random]{rnoise}} \tab random pixel noise
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.

    \tabular{ll}{
      \code{\link[spatstat.explore]{density.psp}} \tab kernel smoothing of line segments\cr
      \code{\link[spatstat.random]{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link[spatstat.random]{rpoislinetess}} \tab generate tessellation using Poisson line
      process 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link[spatstat.random]{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link[spatstat.random]{rpoispp3}} \tab generate Poisson random points in 3-D \cr
      \code{\link[spatstat.explore]{envelope.pp3}} \tab generate simulation envelopes for
      3-D pattern \cr
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link[spatstat.random]{runifpointx}} \tab generate uniform random points \cr
      \code{\link[spatstat.random]{rpoisppx}} \tab generate Poisson random points 
    }
    

  \bold{Classical exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{clarkevans}} \tab Clark and Evans aggregation index \cr
    \code{\link[spatstat.explore]{fryplot}} \tab Fry plot \cr
    \code{\link[spatstat.explore]{miplot}} \tab Morisita Index plot
  }

  \bold{Smoothing:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{density.ppp}} \tab kernel smoothed density/intensity\cr
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.explore]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.explore]{bw.diggle}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{density.ppp}}\cr
    \code{\link[spatstat.explore]{bw.ppl}} \tab likelihood cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{density.ppp}}\cr
    \code{\link[spatstat.explore]{bw.CvL}} \tab Cronie-Van Lieshout bandwidth selection
                                for density estimation\cr
    \code{\link[spatstat.explore]{bw.scott}} \tab Scott's rule of thumb 
                                for density estimation\cr
    \code{\link[spatstat.explore]{bw.abram.ppp}} \tab Abramson's rule for adaptive bandwidths\cr
    \code{\link[spatstat.explore]{bw.relrisk}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{relrisk}} \cr
    \code{\link[spatstat.explore]{bw.smoothppp}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{Smooth.ppp}} \cr
    \code{\link[spatstat.explore]{bw.frac}} \tab bandwidth selection using window geometry\cr
    \code{\link[spatstat.explore]{bw.stoyan}} \tab Stoyan's rule of thumb for bandwidth
                               for \code{\link[spatstat.explore]{pcf}}
  }

  \bold{Modern exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{clusterset}} \tab Allard-Fraley feature detection  \cr
    \code{\link[spatstat.explore]{nnclean}} \tab Byers-Raftery feature detection  \cr
    \code{\link[spatstat.explore]{sharpen.ppp}} \tab Choi-Hall data sharpening \cr
    \code{\link[spatstat.explore]{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.explore]{rho2hat}} \tab Kernel estimate of effect of two covariates\cr
    \code{\link[spatstat.explore]{spatialcdf}} \tab Spatial cumulative distribution function\cr
    \code{\link[spatstat.explore]{roc}} \tab Receiver operating characteristic curve \cr
    \code{\link[spatstat.explore]{sdr}} \tab Sufficient Data Reduction \cr
    \code{\link[spatstat.explore]{thresholdSelect}} \tab optimal
    thresholding of a predictor 
  }

  \bold{Summary statistics for a point pattern:}

  \tabular{ll}{
    \code{\link[spatstat.explore]{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.explore]{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link[spatstat.explore]{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link[spatstat.explore]{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link[spatstat.explore]{Lest}} \tab Besag \eqn{L}-function\cr
    \code{\link[spatstat.explore]{Tstat}} \tab Third order \eqn{T}-function \cr
    \code{\link[spatstat.explore]{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link[spatstat.explore]{pcf}} \tab 	pair correlation function \cr
    \code{\link[spatstat.explore]{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Linhom}} \tab \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{pcfinhom}} \tab pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.explore]{Finhom}} \tab \eqn{F} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Ginhom}} \tab \eqn{G} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Jinhom}} \tab \eqn{J} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localL}} \tab Getis-Franklin neighbourhood density function\cr
    \code{\link[spatstat.explore]{localK}} \tab neighbourhood K-function\cr
    \code{\link[spatstat.explore]{localpcf}} \tab local pair correlation function\cr
    \code{\link[spatstat.explore]{localKinhom}} \tab local \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localLinhom}} \tab local \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localpcfinhom}} \tab local pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.explore]{Ksector}} \tab Directional \eqn{K}-function\cr
    \code{\link[spatstat.explore]{Kscaled}} \tab locally scaled \eqn{K}-function \cr
    \code{\link[spatstat.explore]{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link[spatstat.explore]{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link[spatstat.explore]{envelope}} \tab simulation envelopes for a summary
    function \cr
    \code{\link[spatstat.explore]{varblock}} \tab variances and confidence intervals\cr
                   \tab for a summary function \cr
    \code{\link[spatstat.explore]{lohboot}} \tab bootstrap for a summary function 
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.explore]{plot.fv}} \tab plot a summary function\cr
    \code{\link[spatstat.explore]{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link[spatstat.explore]{harmonise.fv}} \tab make functions compatible \cr
    \code{\link[spatstat.explore]{eval.fasp}} \tab evaluate any expression involving
    an array of functions\cr
    \code{\link[spatstat.explore]{with.fv}} \tab evaluate an expression for a 
    summary function\cr
    \code{\link[spatstat.explore]{Smooth.fv}} \tab apply smoothing to a summary function\cr
    \code{\link[spatstat.explore]{deriv.fv}} \tab calculate derivative of a summary function\cr
    \code{\link[spatstat.explore]{pool.fv}} \tab pool several estimates of a summary function\cr
    \code{\link[spatstat.explore]{density.ppp}} \tab kernel smoothed density\cr
    \code{\link[spatstat.explore]{densityHeat.ppp}} \tab diffusion kernel smoothed density\cr
    \code{\link[spatstat.explore]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.explore]{sharpen.ppp}} \tab data sharpening  \cr
    \code{\link[spatstat.random]{rknn}} \tab theoretical distribution of nearest
    neighbour distance
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} such that \code{marks(X)} is a factor. 
  \tabular{ll}{
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimation of relative risk  \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan test of elevated risk  \cr
    \code{\link[spatstat.explore]{Gcross},\link[spatstat.explore]{Gdot},\link[spatstat.explore]{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{G[i,j], G[i.]} \cr
    \code{\link[spatstat.explore]{Kcross},\link[spatstat.explore]{Kdot}, \link[spatstat.explore]{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{K[i,j], K[i.]} \cr
    \code{\link[spatstat.explore]{Lcross},\link[spatstat.explore]{Ldot}} \tab
    multitype \eqn{L}-functions 
    \eqn{L_{ij}, L_{i\bullet}}{L[i,j], L[i.]} \cr
    \code{\link[spatstat.explore]{Jcross},\link[spatstat.explore]{Jdot},\link[spatstat.explore]{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{J[i,j],J[i.]} \cr
    \code{\link[spatstat.explore]{pcfcross}} \tab
    multitype pair correlation function \eqn{g_{ij}}{g[i,j]} \cr
    \code{\link[spatstat.explore]{pcfdot}} \tab
    multitype pair correlation function \eqn{g_{i\bullet}}{g[i.]} \cr
    \code{\link[spatstat.explore]{pcfmulti}} \tab
    general pair correlation function \cr
    \code{\link[spatstat.explore]{markconnect}} \tab
    marked connection function \eqn{p_{ij}}{p[i,j]} \cr
    \code{\link[spatstat.explore]{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link[spatstat.explore]{Iest}} \tab  multitype \eqn{I}-function\cr
    \code{\link[spatstat.explore]{Kcross.inhom},\link[spatstat.explore]{Kdot.inhom}} \tab
    inhomogeneous counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.explore]{Lcross.inhom},\link[spatstat.explore]{Ldot.inhom}} \tab
    inhomogeneous counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.explore]{pcfcross.inhom},\link[spatstat.explore]{pcfdot.inhom}} \tab
    inhomogeneous counterparts of \code{pcfcross}, \code{pcfdot} \cr
    \code{\link[spatstat.explore]{localKcross},\link[spatstat.explore]{localKdot}} \tab
    local counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.explore]{localLcross},\link[spatstat.explore]{localLdot}} \tab
    local counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.explore]{localKcross.inhom},\link[spatstat.explore]{localLcross.inhom}} \tab
    local counterparts of \code{Kcross.inhom}, \code{Lcross.inhom} 
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type. For numeric marks, there are the
  following functions:
  \tabular{ll}{
    \code{\link[spatstat.explore]{markmean}} \tab smoothed local average of marks \cr
    \code{\link[spatstat.explore]{markvar}} \tab smoothed local variance of marks \cr
    \code{\link[spatstat.explore]{markcorr}} \tab mark correlation function \cr
    \code{\link[spatstat.explore]{markcrosscorr}} \tab mark cross-correlation function \cr
    \code{\link[spatstat.explore]{markvario}} \tab mark variogram \cr
    \code{\link[spatstat.explore]{markmarkscatter}} \tab mark-mark scatterplot \cr
    \code{\link[spatstat.explore]{Kmark}} \tab mark-weighted \eqn{K} function \cr
    \code{\link[spatstat.explore]{Emark}} \tab mark independence diagnostic \eqn{E(r)} \cr
    \code{\link[spatstat.explore]{Vmark}} \tab mark independence diagnostic \eqn{V(r)} \cr
    \code{\link[spatstat.explore]{nnmean}} \tab nearest neighbour mean index \cr
    \code{\link[spatstat.explore]{nnvario}} \tab nearest neighbour mark variance index 
  }
  For marks of any type, there are the following:
  \tabular{ll}{
    \code{\link[spatstat.explore]{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link[spatstat.explore]{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link[spatstat.explore]{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link[spatstat.geom]{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{marktable}} \tab tabulate the marks of neighbours
    in a point pattern 
  }

  \bold{Summary statistics for a three-dimensional point pattern:}

  These are for 3-dimensional point pattern objects (class \code{pp3}).

  \tabular{ll}{
    \code{\link[spatstat.explore]{F3est}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.explore]{G3est}} \tab nearest neighbour function \eqn{G} \cr
    \code{\link[spatstat.explore]{K3est}} \tab \eqn{K}-function \cr
    \code{\link[spatstat.explore]{pcf3est}} \tab pair correlation function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.explore]{envelope.pp3}} \tab simulation envelopes 
  }

  \bold{Summary statistics for random sets:}
  
  These work for point patterns (class \code{ppp}),
  line segment patterns (class \code{psp})
  or windows (class \code{owin}).
  
  \tabular{ll}{
    \code{\link[spatstat.explore]{Hest}} \tab spherical contact distribution \eqn{H} \cr
    \code{\link[spatstat.explore]{Gfox}} \tab Foxall \eqn{G}-function \cr
    \code{\link[spatstat.explore]{Jfox}} \tab Foxall \eqn{J}-function
  }

 \bold{Model fitting}

 Functions for fitting point process models
 are now contained in the \pkg{spatstat.model} package.
  

  \bold{Simulation}

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}
  Functions for random generation are now contained in the \pkg{spatstat.random} package.

  See also \code{\link[spatstat.explore]{varblock}} for estimating the variance
  of a summary statistic by block resampling, and
  \code{\link[spatstat.explore]{lohboot}} for another bootstrap technique.

  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.
   
   Methods for simulating a fitted model are now
   contained in the \pkg{spatstat.model} package.

   \bold{Other random patterns:}
   Functions for random generation are now contained in the
   \pkg{spatstat.random} package.

   \bold{Simulation-based inference}

   \tabular{ll}{
    \code{\link[spatstat.explore]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.explore]{bits.envelope}} \tab critical envelope for balanced
    two-stage Monte Carlo test \cr
    \code{\link[spatstat.model]{qqplot.ppm}} \tab diagnostic plot for interpoint
    interaction \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.explore]{studpermu.test}} \tab studentised permutation test\cr
    \code{\link[spatstat.explore]{segregation.test}} \tab test of segregation of types 
  }

  \bold{Hypothesis tests:}

  \tabular{ll}{
    \code{\link[spatstat.explore]{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link[spatstat.explore]{clarkevans.test}} \tab Clark and Evans test \cr
    \code{\link[spatstat.explore]{cdf.test}} \tab Spatial distribution goodness-of-fit test\cr
    \code{\link[spatstat.explore]{berman.test}} \tab Berman's goodness-of-fit tests\cr
    \code{\link[spatstat.explore]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.explore]{dclf.test}} \tab Diggle-Cressie-Loosmore-Ford test \cr
    \code{\link[spatstat.explore]{mad.test}} \tab Mean Absolute Deviation test \cr
    \code{\link[spatstat.model]{anova.ppm}} \tab Analysis of Deviance for
    point process models 
  }

  \bold{More recently-developed tests:}
  
  \tabular{ll}{
    \code{\link[spatstat.explore]{dg.test}} \tab Dao-Genton test \cr
    \code{\link[spatstat.explore]{bits.test}} \tab Balanced independent two-stage test \cr
    \code{\link[spatstat.explore]{dclf.progress}} \tab Progress plot for DCLF test \cr
    \code{\link[spatstat.explore]{mad.progress}} \tab Progress plot for MAD test \cr
  }

  \bold{Model diagnostics:}

  Classical measures of model sensitivity such as leverage and
  influence, and classical model diagnostic tools such as 
  residuals, partial residuals, and effect estimates,
  have been adapted to point process models. 
  These capabilities are now provided in the \pkg{spatstat.model} package.
  

  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning  
  }
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Julian Gilbey, 
  Yongtao Guan,
  Ute Hahn,
  Kassel Hingee,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  For comments, corrections, bug alerts and suggestions, we thank
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

