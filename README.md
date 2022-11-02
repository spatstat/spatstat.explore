# spatstat.explore

## Exploratory/nonparametric data analysis for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.explore)](https://CRAN.R-project.org/package=spatstat.explore) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.explore)](https://github.com/spatstat/spatstat.explore)

The original `spatstat` package has been split into
several sub-packages (See [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.explore` is one of the
sub-packages. It contains the main user-level functions that perform
**exploratory** and **nonparametric** statistical analysis of spatial data,
with the exception of data on linear networks.

Most of the functionality is for spatial point patterns in two dimensions.
There is a very modest amount of functionality for 3D and higher dimensional patterns
and space-time patterns.

### Overview 

`spatstat.explore` supports

- data manipulation and exploratory graphics
- exploratory analysis 
- smoothing
- cluster detection
- nonparametric estimation 
- hypothesis tests (simulation-based and nonparametric)

### Detailed contents

For a full list of functions, see the help file for `spatstat.explore-package`.

#### Exploratory analysis 

- Clark-Evans index, Hopkins-Skellam index
- quadrat counting estimates of intensity, quadrat counting test
- Fry plot
- Morisita plot
- scan statistic
- cluster detection (Allard-Fraley cluster set, Byers-Raftery cleaning)

#### Nonparametric estimation

- kernel estimation of intensity of a point pattern
- kernel smoothing of mark values attached to point locations
- kernel estimation of relative risk
- kernel smoothing of a line segment pattern
- bandwidth selection
- nonparametric estimation of intensity as a function of a covariate
- ROC curve, AUC
- summary functions (K-function, pair correlation function,
empty space function, nearest neighbour distance function, J-function, etc)
and multi-type versions of these functions
- mark correlation function, mark independence diagnostoc
- local summary functions (LISA)
- simulation envelopes of summary functions
- manipulation of summary functions (plot, evaluate, differentiate, smooth etc)
- spatial bootstrap

#### Formal inference

- hypothesis tests (quadrat test, Clark-Evans test, Berman test, Diggle-Cressie-Loosmore-Ford test, scan test, studentised permutation test, segregation test, envelope tests, Dao-Genton test, balanced independent two-stage test)

#### Data manipulation

- image blurring
- Choi-Hall data sharpening of point locations
- transects of an image along a line or curve
- programming tools

