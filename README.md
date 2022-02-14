<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cort)](https://CRAN.R-project.org/package=cort)
[![CRAN number of downloads](https://cranlogs.r-pkg.org/badges/grand-total/cort)](https://cranlogs.r-pkg.org/badges/grand-total/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
[![DOI](https://zenodo.org/badge/247063359.svg)](https://zenodo.org/badge/latestdoi/247063359)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02653/status.svg)](https://doi.org/10.21105/joss.02653)
[![tic](https://github.com/lrnv/cort/workflows/tic/badge.svg?branch=master)](https://github.com/lrnv/cort/actions)
<!-- badges: end -->

The `cort` package provides S4 classes and methods to fit several copula models: 

* The classic empirical checkerboard copula and the empirical checkerboard copula with known margins, see Cuberos, Masiello and Maume-Deschamps (2019) are proposed. These two models allow to fit copulas in high dimension with a small number of observations, and they are always proper copulas. Some flexibility is added via a possibility to differentiate the checkerboard parameter by dimension. 

* The last model consist of the implementation of the Copula Recursive Tree algorithm, aka. CORT, including the localised dimension reduction, which fits a copula by recursive splitting of the copula domain, see Laverny, Maume-Deschamps, Masiello and Rullière (2020).

* We finally provide an efficient way of mixing copulas, allowing to bag the algorithm into a forest, and a generic way of measuring d-dimensional boxes with a given copula.

## Installation

`cort` is Now on [CRAN](https://CRAN.R-project.org)! You can install the stable version with:

``` r
install.packages("cort")
```

The upstream development version can also be installed with :

``` r
devtools::install_github("lrnv/cort")
```

Note that the installation from github will require the system to have a compiler: 

- Windows: Rtools
- macOS: Xcode CLI
- Linux: r-base-dev (debian)


The vignettes are quite expressive. They give a clear overview of what can be done with this package, how it is coded and why it is useful. Please read them for more details. 

## How to report bugs and get support

To report a bug, feel free to open an issue on the github repository. Support can also be provided through the same chanel if you need it.

## How to contribute

Every contribution is welcome, on the form of pull requests on the github repository. For large modifications, please open an issue for discussions firsts. Concerning the naming convention, the CamelCase functions usually designate classes and constructors of these classes, and all other methods are in snake_case.


## Citation

If you use this work, you may cite the following references. To refer to the theory of the CORT estimator, you may cite : 

```bib
@article{laverny2021dependence,
  title = {Dependence structure estimation using Copula Recursive Trees},
  journal = {Journal of Multivariate Analysis},
  volume = {185},
  pages = {104776},
  year = {2021},
  issn = {0047-259X},
  doi = {10.1016/j.jmva.2021.104776},
  author = {Oskar Laverny and Esterina Masiello and Véronique Maume-Deschamps and Didier Rullière}
}
```

To refer to the package itself, you may cite: 

```bib
@article{laverny2020empirical,
  doi = {10.21105/joss.02653},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {56},
  pages = {2653},
  author = {Oskar Laverny},
  title = {Empirical and non-parametric copula models with the {cort R} package},
  journal = {Journal of Open Source Software}
}
```
