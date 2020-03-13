#' @import methods
#' @import copula
#' @importFrom stats runif
#' @importFrom utils capture.output
NULL



#' Copula volume on hyper-boxes
#'
#' u must be piecewise smaller than v, otherwise the function will return an error.
#'
#' A method is currently implemented for the main virtual class 'Copula', but it assumes
#' that a pCopula method is avaliable for the given copula.
#'
#' This function calculates the measure of the copula according to the algorythme proposed by :
#' Umberto Cherubini & Silvia Romagnoli (2009) Computing the
#' Volume of n-Dimensional Copulas, Applied Mathematical Finance, 16:4, 307-314, DOI:
#'   10.1080/13504860802597311 link : \url{http://dx.doi.org/10.1080/13504860802597311}
#'
#'
#' @param u numeric matrix : minimum point of the hyper-rectangles, one row per observation.
#' @param v numeric matrix : maximum point of the hyper-rectangle, one row per observation.
#' @param copula the copula to calcule it's measure on [u,v]
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the measure of the copula.
#' @exportMethod vCopula
#' @name vCopula
#' @rdname vCopula-methods
#'
#' @examples
#' library(empCop)
#' # For a simple one-dimentional input :
#' cop = copula::archmCopula('Clayton',0.7,3)
#' vCopula(rep(0,3),rep(1,3),cop)
#' # the function is vectorised :
#' v=matrix(seq(0,1,length.out=12),ncol=3)
#' u=matrix(rep(0,12),ncol=3)
#' vCopula(u,v,cop)
setGeneric("vCopula", function(u, v, copula, ...) {

  # taken from the generic of pCopula, does mainly the same...

  if (!is.matrix(u))
    u <- rbind(u, deparse.level = 0L)
  ## here as well, 'outside' and 'on-boundary' are equivalent:
  u[] <- pmax(0, pmin(1, u))

  if (!is.matrix(v))
    v <- rbind(v, deparse.level = 0L)
  ## here as well, 'outside' and 'on-boundary' are equivalent:
  v[] <- pmax(0, pmin(1, v))

  standardGeneric("vCopula")
})
