#' @include generics.R
NULL

.Box = setClass(Class = "Box",
                slots = c(a = "numeric",
                          b="numeric",
                          dim="numeric",
                          volume="numeric",
                          middle_point="numeric"),
                validity = function(object) {
                   errors <- c()
                   if(length(object@a) != object@dim){
                     errors <- c(errors,"you should provide the right dimension at initialisation.")
                   }
                   if(length(object@a) != length(object@b)){
                     errors <- c(errors, "you should provide vectors a and b of same length")
                   }
                   if (!all((0 <= object@a)*(object@a <= object@b)*(object@b <= 1))) {
                     errors <- c(errors, "You should provide a and b s.t 0 <= a <= b <= 1")
                   }
                   if(any(object@middle_point != (object@a + (object@b-object@a)/2))){
                     errors <- c(errors, "The middle point is not ok")
                   }
                   if(!(object@volume == prod(object@b-object@a))){
                     errors <- c(errors,"volume is not OK")
                   }
                   if (length(errors) == 0)
                     TRUE else errors
                 })

Box = function(a,b) {
  if((all(a == b) * all(a == 0)) || (any(a == b))){
    return(ZeroBox(length(a)))
  }
  return(.Box(a = a,
              b = b,
              dim = length(a),
              volume = prod(b-a),
              middle_point = a+(b-a)/2))
}

ZeroBox = function(dim){
  return(.Box(a = rep(0,dim),
              b = rep(0,dim),
              dim = dim,
              volume = 0,
              middle_point = rep(0,dim)))
}

setMethod(f = "show",
          signature = c(object = "Box"),
          definition = function(object){
  cat("Box:[",object@a,"]-> [",object@b,"]")
})

setMethod(f = "intersect",
          signature = c(object = "Box",b="Box"),
          definition = function(object,b){
            if(dim(b) != object@dim){
              stop("boxes must have same dimension")
            }
            if(all((object@a < b@b) * (object@b > b@a))){
              return(Box(pmax(object@a,b@a),pmin(object@b,b@b)))
            } else {
              return(ZeroBox(object@dim))
            }
          })

setMethod(f = "simu_unif",
          signature = c(object = "Box",n="numeric"),
          definition = function(object,n){
            if(n==0){
              return(matrix(0,nrow=0,ncol=object@dim))
            }
            return(t(object@a + matrix(runif(n*object@dim),nrow=object@dim) * (object@b - object@a)))
          })

setMethod(f = "contains",
          signature = c(object = "Box",u="matrix"),
          definition = function(object,u,type="loose"){
            if(type == "strict"){
              return(apply((t(u) > object@a)*(t(u)<object@b),2,all))
            }
            if(type == "loose"){
              return(apply((t(u) >= object@a)*(t(u)<=object@b),2,all))
            }
            if(type == "rllc"){
              return(apply((t(u) >= object@a)*(t(u)<object@b),2,all))
            }
            stop("type should be one of loose, strict of rllc")
          })

setMethod(f = "measure_in",
          signature = c(object = "Box",u="matrix"),
          definition = function(object,u){
            if ((length(u)%%object@dim) != 0) stop("you must provide the right number of dimensions")
            apply(pmax(pmin(t(u),object@b) - object@a,0),2,prod)/object@volume
          })

setMethod(f = "project",
          signature = c(object = "Box",dimensions = "numeric"),
          definition = function(object,dimensions){
            if(!all(dimensions %in% 1:object@dim)){
              stop("you must provide dimensions that re inside the Ok dimensiosn of the box")
            }
            return(Box(object@a[dimensions],object@b[dimensions]))
          })

setMethod(f = "split",
          signature = c(object = "Box"),
          definition = function(object,breakpoint,breakpoint_dim){

            if(length(breakpoint) != length(breakpoint_dim)){
              stop("length of breakpoint should be equal to length of breakpoint_dim")
            }
            if(length(breakpoint) > object@dim){
              stop("length of breakpoint should be smaller than the dimension of the box")
            }
            if(any(!(breakpoint_dim %in% (1:object@dim)))){
              stop("breakpoint dimensions should be inside the box dimensions.")
            }
            if(length(unique(breakpoint_dim)) != length(breakpoint_dim)){
              stop("breakpoint dimensions should be unique")
            }

            projection = project(object,breakpoint_dim)
            if(!contains(projection,breakpoint,type="loose")){stop("The box should loosely contain the breakpoint")}
            if(!contains(projection,breakpoint,type="strict")){
              cat("be carrefull, we are removing a dimensions because the breakpoint touches the boundary")
              dims_to_split = (breakpoint != object@a[breakpoint_dim]) * (breakpoint != object@b[breakpoint_dim])
              breakpoint = breakpoint[dims_to_split]
              breakpoint_dim = breakpoint_dim[dims_to_split]
            } else {
              dims_to_split = 1:length(breakpoint_dim)
            }
            d_split=length(breakpoint_dim)

            1:(2^d_split) %>%
              purrr::map(number2binary,d_split) %>%
              purrr::map(function(bin){
                new_a = object@a
                new_b = object@b
                new_a[breakpoint_dim] = breakpoint^bin * object@a[breakpoint_dim] ^ (1-bin)
                new_b[breakpoint_dim] = breakpoint^(1-bin) * object@b[breakpoint_dim] ^ bin
                return(Box(a = new_a, b = new_b))
              }) %>%
              return
          })
















