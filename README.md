<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cort)](https://CRAN.R-project.org/package=cort)
[![Travis build status](https://img.shields.io/travis/com/lrnv/cort/master?logo=travis&style=flat-square&label=Linux)](https://travis-ci.com/lrnv/cort)
[![Github Build status](https://img.shields.io/github/workflow/status/lrnv/cort/R%20CMD%20Check%20via%20%7Btic%7D?logo=github&label=Github%20build&style=flat-square)](https://github.com/lrnv/cort/actions)
[![AppVeyor build status](https://img.shields.io/appveyor/ci/lrnv/cort?label=Windows&logo=appveyor&style=flat-square)](https://ci.appveyor.com/project/lrnv/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
<!-- badges: end -->

The `cort` package provides S4 classes and methods to fit several copula models: 

* The classic empirical checkerboard copula and the empirical checkerboard copula with known margins, see Cuberos, Masiello and Maume-Deschamps (2019) <doi:10.1080/03610926.2019.1586936> are proposed. These two models allow to fit copulas in high dimension with a small number of observations, and they are always proper copulas. Some flexibility is added via a possibility to differentiate the checkerboard parameter by dimension. 

* The last model consist of the implementation of the Copula Recursive Tree algorithm, aka. CORT, including the localised dimension reduction, which fits a copula by recursive splitting of the copula domain. 

* We finaly provide an efficient way of mixing copulas, allowing to bag the algorithm into a forest, and a generic way of measuring d-dimensional boxes with a given copula.

## Installation

`cort` is Now on [CRAN](https://CRAN.R-project.org)! You can install the sable version with:

``` r
install.packages("cort")
```

The upstream developpement version can also be installed with :

``` r
devtools::install_github("lrnv/cort")
```


The vignettes are quite expressive. They give a clear overview of what can be done with this package, how it is coded and why it is useful. Please read them for more details. 

