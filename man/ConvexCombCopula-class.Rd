% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvexCombCopula.R
\name{ConvexCombCopula-Class}
\alias{ConvexCombCopula-Class}
\alias{ConvexCombCopula}
\title{Convex Combination of copulas.}
\usage{
ConvexCombCopula(copulas, alpha = rep(1, length(copulas)))
}
\arguments{
\item{copulas}{a list of copulas of same dimension}

\item{alpha}{a vector of (positive) weights}
}
\value{
An instance of the \code{ConvexCombCopula} S4 class. The object represent the copula that results from a convex combinaison of other copulas, and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
}
\description{
ConvexCombCopula class
}
\details{
The ConvexCombcopula class is used to build convex combinations of copulas,
with given positives weights. The rCopula and pCopula functions works for
those copulas, assuming they work for the given copulas that we combined
in a convex way.

See the corresponding vignette for more details about the implementation.
}
\examples{
dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
copulas <- list(
  cbCopula(dataset[,2:3],m=10),
  cbCopula(dataset[,2:3],m=5)
)
alpha <- c(1,4)
(cop <- ConvexCombCopula(copulas,alpha))
}
