#' @import methods
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @include utils.R utils-pipe.R
NULL

#' Intersection of 2 boxes
#'
#' Gives the intersection of two boxes
#'
#' @param object an object of class `Box`
#' @param b an other object of class `Box`
#'
#' @return A box object representing the intersection of those two boxes.
#' @exportMethod intersect
#' @name intersect
#' @rdname intersect-methods
#'
#' @examples
#' library(cort)
#' intersect(Box(rep(0,2),rep(1/2,2)),Box(rep(1/4,2),rep(1/2,2)))
setGeneric("intersect",function(object,b){
  standardGeneric("intersect")
})

#' Contains method
#'
#' Answer the question : "Does the box contains those points ? The `type` parameter is used to say if you want the test to be strict on both side, loose on both side, or right-limited and left-continuous via the `rllc` possibility.
#'
#' @param object an object of class `Box`
#' @param u a matrix object, representing points of the right dimension in each row
#' @param type a string, one of c("strict","loose","rllc") (default : loose)
#'
#' @return A vector of boolean values, one for each row of the data `u`
#' @exportMethod contains
#' @name contains
#' @rdname contains-methods
#'
#' @examples
#' library(cort)
#' contains(Box(rep(0,2),rep(1/2,2)),rep(2,2))
setGeneric("contains",signature = c("object","u"),
           function(object,u,type="loose"){
  u <- normalise_data(u,object@dim)
  standardGeneric("contains")
})


#' Measure of a point inside a box
#'
#' compute the formula lambda(0,u inter b) / lambda (b) for a box b and a point u
#'
#' @param object an object of class `Box`
#' @param u a matrix object, representing points of the right dimension in each row
#'
#' @return A vector of values between 0 and 1, representing the measure of each point inside the box.
#' @exportMethod measure_in
#' @name measure_in
#' @rdname measure_in-methods
#'
#' @examples
#' library(cort)
#' measure_in(Box(rep(0,2),rep(1/2,2)),matrix(runif(100),nrow=50,ncol=2))
setGeneric("measure_in",function(object,u){
  u <- normalise_data(u,object@dim)
  standardGeneric("measure_in")
})

#' Simulate uniformly from a box
#'
#' simulate random points uniformly inside a box.
#'
#' @param object an object of class `Box`
#' @param n the number of points to simulate
#'
#' @return a matrix with each random point in a row.
#' @exportMethod simu_unif
#' @name simu_unif
#' @rdname simu_unif-methods
#'
#' @examples
#' library(cort)
#' simu_unif(Box(rep(0,2),rep(1/2,2)),10)
setGeneric("simu_unif",function(object,n){
  standardGeneric("simu_unif")
})


#' project a box onto smaller dimensions
#'
#' Given new dimensions, a box can be projected on these dimensions.
#'
#' @param object an object of class `Box`
#' @param dimensions the dimensions to project on
#'
#' @return a matrix with each random point in a row.
#' @exportMethod project
#' @name project
#' @rdname project-methods
#'
#' @examples
#' library(cort)
#' project(Box(rep(0,2),rep(1/2,2)),1)
setGeneric("project",function(object,dimensions){
  standardGeneric("project")
})

#' Check if a box is splittable
#'
#' To be splitted, a box needs to have at least min_node_size points and 2 splitting dimensions.
#'
#' @param object an object of class `Box`
#' @return a boolean telling you if the box is splittable or not
#' @exportMethod is_splittable
#' @name is_splittable
#' @rdname is_splittable-methods
#'
#' @examples
#' library(cort)
#' b = Box(rep(0,2),rep(1,2))
#' b = WeightedBox(b,matrix(0.5,nrow=2,ncol=2),1,c(1,2))
#' is_splittable(b)
setGeneric("is_splittable",function(object){
             standardGeneric("is_splittable")
           })



#' Split a box onto a given breakpoint
#'
#' This functions splits a box onto given dimensions and given breakpoint
#'
#' @param object an object of class `Box`
#' @param breakpoint a point to split the box on, should be inside the box
#' @param breakpoint_dim the dimensions where you want the splitting to happened.
#' @param ... other parameters that will be passed to methods.
#' @return a matrix with each random point in a row.
#' @exportMethod split
#' @name split
#' @rdname split-methods
#'
#' @examples
#' library(cort)
#' split(Box(rep(0,2),rep(1/2,2)),rep(1/4,2),c(1,2))
setGeneric("split",signature = c("object"),
           function(object,...){
  standardGeneric("split")
})

#' Fit the tree
#'
#' Fits the tree
#'
#' @param object the tree to be fitted
#'
#' @return nothing
#' @exportMethod fit
#' @name fit
#' @rdname fit-methods
setGeneric("fit",function(object){standardGeneric("fit")})





