#' @import methods
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @include utils.R utils-pipe.R
NULL

setGeneric("intersect",function(object,b){
  standardGeneric("intersect")
})

setGeneric("contains",signature = c("object","u"),
           function(object,u,type="loose"){
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

setGeneric("is_splittable",function(object){
             standardGeneric("is_splittable")
           })

setGeneric("split",signature = c("object"),
           function(object,...){
  standardGeneric("split")
})

setGeneric("fit",function(object){standardGeneric("fit")})





