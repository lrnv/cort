#' @import methods
#' @importFrom Rdpack reprompt
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @include utils.R utils-pipe.R
NULL

########## Non-exported generics :


setGeneric("intersect",function(object,b){standardGeneric("intersect")})

setGeneric("is_splittable",function(object){standardGeneric("is_splittable")})

setGeneric("split",signature = c("object"),function(object,...){standardGeneric("split")})

setGeneric("fit",function(object){standardGeneric("fit")})

setGeneric("contains",signature = c("object","u"),function(object,u,type="loose"){
  u <- normalise_data(u,object@dim)
  standardGeneric("contains")
})

setGeneric("measure_in",function(object,u){
  u <- normalise_data(u,object@dim)
  standardGeneric("measure_in")
})

setGeneric("simu_unif",function(object,n){
  standardGeneric("simu_unif")
})

setGeneric("project",function(object,dimensions){
  standardGeneric("project")
})

##########" Exported generics :

#' Copula volume on hyper-boxes
#'
#' u must be piecewise smaller than v, otherwise the function will return an error.
#'
#' A method is currently implemented for the main virtual class 'Copula', but it assumes
#' that a pCopula method is avaliable for the given copula.
#'
#' This function computes the measure of the copula according to the algorithm proposed by the referenced paper.
#'
#'
#' @param u numeric matrix : minimum point of the hyper-rectangles, one row per observation.
#' @param v numeric matrix : maximum point of the hyper-rectangle, one row per observation.
#' @param copula the copula that we compute the measure on the box (u,v)
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the measure of the copula.
#' @exportMethod vCopula
#' @name vCopula
#' @rdname vCopula-methods
#'
#' @examples
#' library(cort)
#' # The exmeple needs to be re-done.
#'
#' @references
#' \insertRef{cherubini2009}{cort}
#'
setGeneric("vCopula", function(u, v, copula, ...) {

  # taken from the generic of pCopula, does mainly the same...

  u <- normalise_data(u,copula@dim)
  v <- normalise_data(v,copula@dim)
  standardGeneric("vCopula")
})


#' Copula density
#'
#' This function returns the density of a given copula on given observations.
#'
#'
#' @param u numeric matrix : one row per observation
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the density of the copula on each observation
#' @exportMethod dCopula
#' @name dCopula
#' @rdname dCopula-methods
#'
#' @examples
#' library(cort)
#' # The exmeple needs to be re-done.
setGeneric("dCopula", function(u, copula, ...) {
  u <- normalise_data(u,copula@dim)
  standardGeneric("dCopula")
})

#' Copula density
#'
#' This function returns the value of the copula itself on given points.
#'
#'
#' @param u numeric matrix : one row per observation
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the density of the copula on each observation
#' @exportMethod pCopula
#' @name pCopula
#' @rdname pCopula-methods
#'
#' @examples
#' library(cort)
#' # The exmeple needs to be re-done.
setGeneric("pCopula", function(u, copula, ...) {
  u <- normalise_data(u,copula@dim)
  standardGeneric("pCopula")
})

#' @rdname vCopula-methods
#' @aliases vCopula,matrix,matrix,Copula
setMethod("vCopula", signature = c(u = "matrix", v = "matrix"),
          definition = function(u, v, copula) {

            # can handle any copula thant pCopula could handle.
            # shoul be better vetorised...
            # u and v must be numeric, copula must be a copula, and v must be
            # smaller than u
            if (nrow(u) != nrow(v)) {
              stop("u and v must have same shape (same number of row and columns)")
            }
            if((nrow(u) == 0) || (nrow(v) == 0)){
              return(numeric())
            }

            if (any(v < u)) {
              stop("u must be smaller than v !")
            }

            # fastening ? maybe not...
            pCop_method <- selectMethod(pCopula, c("matrix", class(copula)))


            d = dim(copula)
            p <- t(sapply(1:(2^d),function(i){number2binary(i-1,d)}))
            sign <- (-1)^rowSums(p)
            return(sapply(1:nrow(u),function(i){
              if(all(u[i,] == v[i,])){return(0)}
              eval_points <-t(t(p) * as.vector(u[i,]) + t(1-p) * as.vector(v[i,]))
              return(sum(sign * pCopula(eval_points,copula)))
            }))
          })


#' Copula random variables simulation
#'
#' This function simulate random variables from a copula.
#'
#'
#' @param n the number of simulations
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the density of the copula on each observation
#' @exportMethod rCopula
#' @name rCopula
#' @rdname rCopula-methods
#'
#' @examples
#' library(cort)
#' # The exemples needs to be re-done.
setGeneric("rCopula", function(n, copula, ...) standardGeneric("rCopula"))
