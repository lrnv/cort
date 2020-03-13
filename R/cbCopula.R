#' @include generics.R empiricalCopula.R
NULL

############################### Checkerboard copula class #######
.cbCopula = setClass(Class = "cbCopula", contains = "empiricalCopula",
  slots = c(m = "numeric"), validity = function(object) {
   errors <- c()
   if (any(sapply(object@m,function(m){(nrow(object@pseudo_data)%%m) != 0}))) {
     errors <- c(errors, "m should divide the number of row")
   }
   if(length(object@m) != ncol(object@pseudo_data)){
     errors <- c(errors, "lengths of m parameter shoudl be the same as the number of columns of pseudo_data")
   }
   if (length(errors) == 0)
     TRUE else errors
  })
#' cbCopula contructor
#'
#' @param x the data to be used
#' @param m checkerboard parameters
#' @param pseudo Boolean, defaults to `FALSE`. Set to `TRUE` if you are already
#'  providing pseudo datas into the `x` argument.
#'
#' @details The checkerboard copula is a kind of patchwork copula that only uses independant copula as fill-in, only where there are values on the empirical data provided. To create such a copula, you should provide data and checkerboard parameters (depending on the dimention of the data). To read more about thoose models, we refer to the work of Durante & Al, but also on the work of Cuberos & Al
#'
#' @return a cbCopula object
#' @export
cbCopula = function(x, m = rep(nrow(x),ncol(x)), pseudo = FALSE) {
  if (missing(x)) { stop("The argument x must be provided") }
  if(ncol(x) == 0){ stop("you are providing a data.frame equal to NULL") }
  if(nrow(x) == 0){ return(indepCopula(ncol(x))) }
  if (!pseudo) {
    x <- apply(x, 2, rank, na.last = "keep")/(nrow(x) + 1)
  }
  if(length(m) == 1){m = rep(m,ncol(x))}
  if(length(m) != ncol(x)){stop("You should provide m values same lengths as the number of columns in data.")}
  return(.cbCopula(pseudo_data = as.data.frame(x), m = m))
}
setMethod(f = "show",    signature = c(object = "cbCopula"),                definition = function(object)    {
  cat("This is a cbCopula , with : \n", "  dim =", dim(object), "\n   n =",
      nrow(object@pseudo_data), "\n   m =", object@m, "\n")
  cat("The variables names are : ", colnames(object@pseudo_data))
})
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "cbCopula"), definition = function(n, copula) {

  # get info from the copula :
  x = copula@pseudo_data
  x = as.matrix(x)
  rownames(x) <- NULL
  d = ncol(x)
  m = copula@m

  # if n=0, return a 0xdim matrix :
  if (n == 0) { return(matrix(0, nrow = 0, ncol = d)) }

  # Then, let's sample rows coresponding to observations, i.e let's
  # sample thoose boxes with probabilities proportional to the number of
  # observations inside the box.  notes that boxes with probability 0,
  # i.e without observation, were not included here.  This makes the
  # algorythme fast.
  rows <- sample(x = 1:nrow(x), size = n, replace = TRUE)
  seuil_inf <- boxes_from_points(x,m)[rows,,drop=FALSE]
  seuil_sup <- t(t(seuil_inf) + 1/m)

  # Finaly, sample some random uniform, and bound them inside the sampled boxes :
  rng <- matrix(runif(d * n), nrow = n, ncol = d)
  result <- seuil_inf + rng * (seuil_sup - seuil_inf)
  return(result)
})
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "cbCopula"), definition = function(u, copula) {

  # remind that pCopula and dCopula generics already transform inputs
  # into matrices...

  if (ncol(u) != dim(copula)) {
    stop("the input value must be coer??able to a matrix with dim(copula) columns.")
  }

  x = copula@pseudo_data
  x = as.matrix(x)
  rownames(x) <- NULL
  d = ncol(x)
  m = copula@m
  n = nrow(copula@pseudo_data)
  seuil_inf <- boxes_from_points(x,m)

  rez <- sapply(1:nrow(u), function(i) {
    ponderation <- t(apply(seuil_inf, 1, function(y) { u[i, ] - y }))
    ponderation <- pmax(pmin(ponderation, 1/m), 0) ############### Attention probablement une erreur ici dans la vectorisation du m. Cette version est celle avec le m pas vectorisé !!

    sum(apply(ponderation, 1, function(v) { (prod(v) * prod(m)) }))/n
  })

  return(rez)
})
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="cbCopula"),   definition = function(u, copula) {
  stop("Checkerboard copula has no density")
})



