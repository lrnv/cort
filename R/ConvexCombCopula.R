#' @include generics.R
NULL


############################### ConvexComCopula class #######
.ConvexCombCopula = setClass(Class = "ConvexCombCopula", contains = "Copula",
  slots = c(copulas = "list", alpha = "numeric"), validity = function(object) {
    errors <- c()
    if (length(object@copulas) != length(object@alpha)) {
      errors <- c(errors, "the weights parameter alpha must have same length as the copulas list")
    }
    if (!all(sapply(object@copulas, function(x) {
      is(x, "Copula")
    }))) {
      errors <- c(errors, "parameter copulas should contains a list of copulas")
    }
    if (!all(object@alpha >= 0)) {
      errors <- c(errors, "weights should be positives")
    }
    if(abs(sum(object@alpha)-1)>10^(-7)){
      errors <- c(errors,"weights should add up to 1.")
    }
    if (!(length(unique(sapply(object@copulas, dim))) == 1)) {
      errors <- c(errors, "all copulas must have same dimension")
    }
    if (length(errors) == 0)
      TRUE else errors
  })

#' ConvexCombCopula class
#'
#' The ConvexCombcopula class is used to build convex combinations of copulas,
#' with given positives weights. The rCopula and pCopula functions works for
#' thoose copulas, assuming they work for the given copulas that we combined
#' in a convex way.
#'
#' @param copulas a list of copulas of same dimention
#' @param alpha a vector of (positive) weights
#' @name ConvexCombCopula-class
#' @title Convex Combination of copulas.
#'
#' @return a ConvexCombCopula object
#' @export
#'
#' @examples
#' library(empCop)
#' copulas <- list(
#'   copula::archmCopula('gumbel',3,dim=2),
#'   copula::archmCopula('clayton',-1,dim=2)
#' )
#'
#' alpha <- c(1,4)
#'
#' (cop <- ConvexCombCopula(copulas,alpha))
ConvexCombCopula = function(copulas, alpha = rep(1, length(copulas))) {
  if (missing(copulas) || (!is(copulas, "list"))) {
    if(!is(copulas,"Copula")){
      stop("The argument copulas must be provided as a list of copulas")
    } else {
      warning("The copulas argument you provided is a Copula. Returning this copula")
      return(copulas)
    }
  }
  if(length(copulas) == 1){
    warning("You provided is a list of only one copula. Returning this copula")
    return(copulas[[1]])
  }
  .ConvexCombCopula(copulas = copulas, alpha = alpha/sum(alpha))
}

#' @describeIn ConvexCombCopula dimension
#' @param x ConvexCombCopula object
setMethod(f = "dim",     signature = (x = "ConvexCombCopula"),                      definition = function(x)         {
  return(dim(x@copulas[[1]]))
})
setMethod(f = "show",    signature = c(object = "ConvexCombCopula"),                definition = function(object)    {
  cat("This is a ConvexCombCopula , with : \n", "  dim =", dim(object@copulas[[1]]),
      "\n   number of copulas =", length(object@copulas), "\n   alpha =",
      object@alpha, "\n")
  cat("sub-copulas can be accessed trhough the @copulas slot")
})
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "ConvexCombCopula"), definition = function(n, copula) {

  # if n=0, return a 0xdim(copula) matrix :
  if (n == 0) {
    return(matrix(NA, nrow = 0, ncol = dim(copula)))
  }

  # to choose wich copulas will be simulated from, sample
  # 1:length(copulas) with weights equal to alpha, with replacement OFC
  n_cop = length(copula@copulas)
  sampled_copulas <- sample(1:n_cop, size = n, replace = TRUE, prob = copula@alpha)
  tbl <- table(sampled_copulas)

  # then sample from each of thoose copulas the right number of times :
  samples <- mapply(function(x,y){
    rCopula(n = y, copula = copula@copulas[[x]])
  },as.numeric(names(tbl)),as.vector(tbl),SIMPLIFY = FALSE)

  # then rbind all of them and resample (randomly) rows :
  samples <- do.call(rbind, samples)
  samples <- samples[sample(1:nrow(samples), size = nrow(samples),
                            replace = FALSE), ]
  return(samples)
})
setMethod(f = "pCopula", signature = c(u = "matrix", copula = "ConvexCombCopula"),  definition = function(u, copula) {

  # remind that pCopula and dCopula generics already transform inputs into matrices...
  if (ncol(u) != dim(copula)) {
    stop("the input value must be coercable to a matrix with dim(copula) columns.")
  }
  as.vector(
    vapply(copula@copulas,
           pCopula,
           FUN.VALUE=numeric(nrow(u)),
           u=u
    ) %*% copula@alpha
  )

})

