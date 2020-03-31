<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/cort)](https://CRAN.R-project.org/package=cort) -->
[![Travis build status](https://img.shields.io/travis/com/lrnv/cort/master?logo=travis&style=flat-square&label=Linux)](https://travis-ci.com/lrnv/cort)
[![Github Build status](https://img.shields.io/github/workflow/status/lrnv/cort/R%20CMD%20Check%20via%20%7Btic%7D?logo=github&label=Github%20build&style=flat-square)](https://github.com/lrnv/cort/actions)
[![AppVeyor build status](https://img.shields.io/appveyor/ci/lrnv/cort?label=Windows&logo=appveyor&style=flat-square)](https://ci.appveyor.com/project/lrnv/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
<!-- badges: end -->


The `cort` package provides S4 classes to deal with certain non-parametrical copula models : checkerboard constructions, plain or with known margins, convex mixtures of copula models, and the CORT algorithm. This algorithm mimics the classical CART regression estimator to estimate copula, producing a piecewise constant density estimator that is both non-parametric, convergent, and memory-efficient. It includes a localised dimension reduction procedure. Last but not least, we provide a way of bagging copula models that is quite convenient with the CORT algorithm, producing a copula forest. 

## Installation

In the future, you we be able to install the released version of cort from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("cort")
```

but for now, only the upstream version is avaliable. It can be installed with : 

``` r
devtools::install_github("lrnv/cort")
```


The vignettes are quite expressive. They give a clear overview of what can be done with this package, how it is coded and why it is usefull. Please read them for more details. 

