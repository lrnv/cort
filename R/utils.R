#' @include generics.R
NULL

# To include Rcpp, use thoose 2 lines :
# @useDynLib mypkg, .registration = TRUE
# @importFrom Rcpp sourceCpp

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  binary_vector[-(1:(length(binary_vector) - noBits))]
}

# A `sample` function more efficient (cf ?sample,
# exemples)
resample <- function(x, ...) x[sample.int(length(x), ...)]

#' @rdname vCopula-methods
#' @aliases vCopula,matrix,matrix,Copula
setMethod("vCopula", signature = c(u = "matrix", v = "matrix", copula = "Copula"),
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





# box_from_point
# expects a point or several points as a 1xd or nxd matrix and a vector of m's
# DO NOT EXPECT A VECTOR IF THERE IS ONLY ONE POINT
boxes_from_points <- function(points,m){
  t(floor(t(points)*m)/m)
}

















