% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{kendall_func}
\alias{kendall_func}
\alias{kendall_func,Cort-method}
\title{Kendall function of a copula (if it has one)}
\usage{
kendall_func(object, t, ...)

\S4method{kendall_func}{Cort}(object, t, M = 1000)
}
\arguments{
\item{object}{: the tree}

\item{t}{: the value where to compute the kendall function, may be a vector of evaluation values;}

\item{...}{other parameters passed to methods}

\item{M}{the number of simulations}
}
\value{
the quadratic product between the trees
}
\description{
Currently only implemented for Cort models.
Compute the Kendall cdf from the model in a point t
}
\section{Functions}{
\itemize{
\item \code{kendall_func,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
kendall_func(cop,0.5)

}
