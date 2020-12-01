---
title: 'Empirical and non-parametric copula models with the `cort` R package'
tags:
  - R
  - copula
  - statistics
  - ranks
authors:
  - name: Oskar Laverny^[Corresponding author]
    orcid: 0000-0002-7508-999X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Institut Camille Jordan, Université Lyon 1, Lyon, France
   index: 1
 - name: SCOR SE, Paris, France
   index: 2
date: 13 August 2017
bibliography: paper.bib
---

# Summary

The R package `cort` implements object-oriented classes and methods to estimate, simulate and visualize certain types of non-parametric copulas.


Copulas are functions that describe dependence structure of a given dataset, or of a given multivariate random event, without describing the univariate events, the marginals. In statistics, it is sometimes useful to separate the marginal distributions (inflation of the money, mortality of the population) from the dependence structure between them, since estimating everything separately is usually easier. Copulas are broadly used in finance, actuarial science, geostatistics, biostatistics, and many other fields, when dealing with dependence.

Copulas are distributions functions on the unit hypercube that have uniform margins (what we call the 'copula constraints'), and hence this package can be classified in 'density estimation software'. Although the estimation of copulas is a wide-treated subject, most performing estimators available in the literature are based on restricted, parametric estimation: vine copulas [@nagler2016evading] and graphical models [@li2018panda] for example are potential solutions but under restrictive assumptions. Classical density estimators such as kernels or wavelets do not satisfy marginal copula constraints. There also exist several tree-structured piecewise constant density estimators, but they do not always lead to proper copulas when applied on pseudo-observations or true copula samples. The new models that are implemented in this package try to solve these issues.

We note that a lot of tools are available in R for copula modeling  through the excellent package `copula` [@cop1], [@cop2], [@cop3] and [@cop4]. Most of these tools however focus on parametric estimation. We start to bridge the gap by providing some tools for non-parametric estimation.


# Statement of need 

The Copula recursive tree, a.k.a Cort, designed by [@laverny2020dependence] is a flexible, consistent, piecewise linear estimator for a copula [@sklar1959fonctions], leveraging the patchwork copula formalization [@durante2015convergence] and a specific piecewise constant density estimator, the density estimation tree [@ram2011density]. While the patchwork structure imposes the grid, this estimator is data-driven and constructs the grid recursively from the data, minimizing a chosen distance on the copula space. Furthermore, while the addition of the copula constraints makes the available solutions for density estimation unusable, our estimator is only concerned with dependence and guarantees the uniformity of margins. The R package `cort` provides a useful implementation of this model and several potential refinements, allowing for fast computations of Cort trees, and parallel computations of Cort forests.

The main feature implemented in the package is the Cort algorithm, a non-parametric, piecewise constant, copula density estimator. The implementation is recursive and hence quite efficient. The `cort` package is a statistical package that allows to estimate several non-parametric copula models in R. Although the state of the art `copula` package has functions to estimate the empirical copula, we provide a structured set of S4 classes that allows estimation of empirical copulas, checkerboard copulas, Cort copula and bagging of all of these. A specific class exists for bagging Cort models, which implementation runs in parallel, to fasten the computations, using the `future` package [@future]. Most of the underlying machinery and computations are written in `C++`, through the `Rcpp` [@rcpp1], [@rcpp2], [@rcpp3] and [@rcpp4] package.

The `cort` package was designed to be used by statisticians who need a non-parametric view of the dependence structure of a given dataset. It features a rich and extensive API to call the statistic fitting procedures, plotting functions and tools to assess the quality of the fits, while complying with the R standards. Examples datasets are included in the package, and the many vignettes give examples of use cases.

# Acknowledgements

We acknowledge contributions from Véronique Maume-deschamps, Didier Rullière and Esterina Masiello, who all gave meaningful insights about performances of the code. We are thankful to the reviewers for the fruitful review of this article and of the package code.

# References
