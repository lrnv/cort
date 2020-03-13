#' @include generics.R
NULL


############################### Empirical copula class ######
#' Empirical Copula class (virtual mother class)
#'
#' @slot pseudo_data matrix : pseudo_data that the empirical copula is based on.
#'
#'
#' @export
setClass(Class = "empiricalCopula", contains = c("VIRTUAL", "Copula"),
   slots = c(pseudo_data = "data.frame"), validity = function(object) {
     errors <- c()
     if(ncol(object@pseudo_data) == 0){
       errors <- c(errors, "you are providing a data.frame equal to NULL")
     }
     if (prod(apply(object@pseudo_data, 1:2, is.numeric)) != 1) {
       errors <- c(errors, "the data argument must be a numeric data.frame")
     }
     if (prod(object@pseudo_data <= 1) * prod(object@pseudo_data >=
                                              0) == 0) {
       errors <- c(errors, "the pseudo-data should be numeric between 0 and 1 (both included)")
     }
     if (nrow(object@pseudo_data) == 0){
       errors <- c(errors,"pseudo_data data.frame must have at least one row.")
     }
     if (length(errors) == 0)
       TRUE else errors
   })

#' @describeIn empiricalCopula dimension
#' @param x ConvexCombCopula object
setMethod(f = "dim", signature = (x = "empiricalCopula"), definition = function(x) {
  return(ncol(x@pseudo_data))
})

