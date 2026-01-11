# spatstat.explore

## Exploratory/nonparametric data analysis for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.explore)](https://CRAN.R-project.org/package=spatstat.explore) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.explore)](https://github.com/spatstat/spatstat.explore)

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.explore`.
For the latest public release on CRAN, click the green badge above.

 - [Overview of `spatstat.explore`](#overview)
 - [Detailed contents of package](#detailed)
 - [Installing the package](#installing)
 - [Bug reports](#bugreports)
 - [Questions](#questions)
 - [Proposing changes to code](#proposing)
 - [Future development](#future)

___

### <a name="overview"></a> Overview of `spatstat.explore`

The original `spatstat` package has been split into
several sub-packages (See [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.explore` is one of the
sub-packages. It contains the main user-level functions that perform
**exploratory** and **nonparametric** statistical analysis of spatial data,
with the exception of data on linear networks.

Most of the functionality is for spatial point patterns in two dimensions.
There is a very modest amount of functionality for 3D and higher dimensional patterns
and space-time patterns.

`spatstat.explore` supports

- data manipulation and exploratory graphics
- exploratory analysis 
- smoothing
- cluster detection
- nonparametric estimation 
- hypothesis tests (simulation-based and nonparametric)

___

### <a name="detailed"></a> Detailed contents of `spatstat.explore`

For a full list of functions, see the help file for `spatstat.explore-package`.

#### Exploratory analysis 

- Clark-Evans index, Hopkins-Skellam index
- quadrat counting estimates of intensity, quadrat counting test
- Fry plot
- Morisita plot
- scan statistic
- cluster detection (Allard-Fraley cluster set, Byers-Raftery cleaning)

#### Nonparametric estimation of trend

- kernel estimation of intensity of a point pattern
- kernel smoothing of mark values attached to point locations
- kernel estimation of relative risk
- kernel smoothing of a line segment pattern
- bandwidth selection
- spatial CDF 
- nonparametric estimation of intensity as a function of a covariate
- ROC curve, AUC
- Sufficient Data Reduction
- optimal thresholding of a covariate

#### Nonparametric estimation of dependence between points

- summary functions (K-function, pair correlation function,
empty space function, nearest neighbour distance function, J-function, etc)
and multi-type versions of these functions
- mark correlation function, mark independence diagnostic
- local summary functions (LISA)
- simulation envelopes of summary functions
- manipulation of summary functions (plot, evaluate, differentiate, smooth etc)

#### Formal inference

- spatial bootstrap
- asymptotic variance estimates
- hypothesis tests (quadrat test, Clark-Evans test, Berman test, Diggle-Cressie-Loosmore-Ford test, scan test, studentised permutation test, segregation test, envelope tests, Dao-Genton test, balanced independent two-stage test)

#### Data manipulation

- image blurring
- Choi-Hall data sharpening of point locations
- transects of an image along a line or curve
- programming tools

___

### <a name="installing"></a> Installing the package

This repository contains the _development version_ of
`spatstat.explore`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.explore", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.explore`,
type

```R
install.packages("spatstat.explore")
```

___

## <a name="bugreports"></a> Bug reports 

Users are encouraged to report bugs.
If you find a bug in a `spatstat` function,
please identify the sub-package containing that function.
Visit the GitHub repository for the sub-package, 
click the `Issues` tab at the top of the page, 
and press *new issue* to start a new bug report, documentation correction
or feature request.

**Please do not post questions** on the Issues pages,
because they are too clunky for correspondence.

## <a name="questions"></a> Questions about spatstat

For questions about the `spatstat` package family, first check 
the question-and-answer website
[stackoverflow](http://stackoverflow.com/questions/tagged/spatstat)
to see whether your question has already been asked and answered.
If not, you can either post your question at stackoverflow, or
email the authors.

## <a name="proposing"></a> Proposing changes to the code

Feel free to fork `spatstat.explore`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

## <a name="future"></a> Future development

The `spatstat` package family is the result of 30 years of software development
and contains over 200,000 lines of code.
It is still under development,
motivated by the needs of researchers in many fields,
and driven by innovations in statistical science.
We welcome contributions of code, and suggestions
for improvements.
