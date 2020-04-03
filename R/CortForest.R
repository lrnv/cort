.CortForest = setClass(Class = "CortForest", contains = "empiricalCopula",
                 slots = c(p_value_for_dim_red = "numeric",
                           number_max_dim = "numeric",
                           verbose_lvl="numeric",
                           trees = "list",
                           weights = "numeric",
                           indexes = "matrix",
                           pmf = "matrix",
                           norm_matrix = "matrix",
                           oob_pmf = "matrix",
                           oob_kl = "numeric",
                           oob_ise = "numeric"),
                 validity = function(object) {
                   errors <- c()
                   if(!(length(dim(object@data))==2)){
                     errors <- c(errors, "data should be a matrix")
                   }
                   if(!(all(object@data >= 0) & all(1 > object@data))){
                     errors <- c(errors,"data should be a matrix inside 0,1")
                   }
                   if(object@number_max_dim > ncol(object@data)){
                     errors <- c(errors, "the maximum number of splitting dimensions should be smaller than the number of dimensons!")
                   }
                   if (length(errors) == 0)
                     TRUE else errors
                 })

#' CortForest class
#'
#' This class implements the bagging of CORT models, with an oob error minimisation in the weights.
#'
#'
#' @param x The data, must be provided as a matrix with each row as an observation.
#' @param p_value_for_dim_red a p_value for the localised dimension reduction test
#' @param min_node_size The minimum number of observation avaliable in a leaf to initialise a split.
#' @param pseudo_data set to True if you are already providing data on the copula space.
#' @param compte_loo_weights Defaults to FALSE. Allows to use an automatic re-weighting of the trees in the forest, based on leave-one-out considerations.
#' @param number_max_dim The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.
#' @param n_trees Number of trees
#' @param verbose_lvl verbosity level : can be 0 (none) or an integer. bigger the integer bigger the output level.
#'
#' @name CortForest-Class
#' @title Bagged Cort estimates
#' @rdname CortForest-Class
#'
#' @return a CortForest object that can be fitted easily to produce a copula estimate.
#' @export
#'
#' @examples
#' (CortForest(LifeCycleSavings[,1:3],number_max_dim=2,n_trees=2))
CortForest = function(x,
                p_value_for_dim_red=0.75,
                n_trees = 10,
                compte_loo_weights = FALSE,
                min_node_size=1,
                pseudo_data=FALSE,
                number_max_dim=NULL,
                verbose_lvl=2) {

  # coerce the data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL

  if(!pseudo_data){
    data = apply(data,2,rank,ties.method="random")/(nrow(data)+1)
  }

  d = ncol(data)
  n = nrow(data)

  # for each tree, indexes of it's values :
  indexes = matrix(resample(1:n,size=n*n_trees,replace = TRUE),nrow=n_trees,ncol=n) # n_trees, n


  if(verbose_lvl <= 1){
    affichage = ""
  } else {
    affichage = "========================"
  }

  if(verbose_lvl>0){cat(affichage,"Computing trees...\n")}
  trees = furrr::future_map(1:n_trees,function(i){
    return(Cort(data[indexes[i,],],
         p_value_for_dim_red=p_value_for_dim_red,
          min_node_size=min_node_size,
          pseudo_data=TRUE,
          number_max_dim=number_max_dim,
          verbose_lvl=0))},.progress=TRUE)

    if(verbose_lvl>0){cat(affichage,"Computing statistics...\n")}
  # now compute the masked pmf :
  if(verbose_lvl>1){cat("     Computing pmf...\n")}
  pmf = matrix(0,nrow=n_trees,ncol=n)
  is_in = matrix(1:n,nrow=n,ncol=n_trees) # n, n_trees
  for (i in 1:n_trees){
    for (j in 1:n){
      is_in[j,i] = (j %in% indexes[i,])
    }
    pmf[i,] = dCopula(data,trees[[i]])
  }

  if(verbose_lvl>1){cat("     Computing norm matrix...\n")}
  norm_matrix = normMatrix(as = purrr::map(trees,~.x@a),
                           bs = purrr::map(trees,~.x@b),
                           kernels = purrr::map(trees,~.x@p/.x@vols))

  # we now need the weights. Let's not fit them yet.
  weights = rep(1:n_trees,n_trees)

  # OUT OF BAG STATS
  if(verbose_lvl>1){cat("     Computing oob stats...\n")}
  oob_pmf = matrix(0,nrow=n,ncol=n_trees-1)
  oob_kl = numeric(n_trees-1)
  oob_ise = numeric(n_trees - 1)

  for (j in 2:n_trees){
    weights = rep(1/j,j) # weights computed from only theese trees

    oob_wts = weights
    dim(oob_wts) = c(1,j)
    oob_wts = t(oob_wts[rep(1,n),]*(1-is_in[,1:j]))
    oob_pmf[,j-1] = colSums(pmf[1:j,]*oob_wts)/colSums(oob_wts)
    oob_kl[j-1] = -mean(log(oob_pmf[!is.na(oob_pmf[,j-1]),j-1]))
    oob_ise[j-1] = weights %*% norm_matrix[1:j,1:j] %*% weights
  }

  if(verbose_lvl>0){cat(affichage,"Done !\n")}
  # Output the forest with all it's stats :
  return(.CortForest(
    data = data,
    p_value_for_dim_red = p_value_for_dim_red,
    number_max_dim = max(number_max_dim,d),
    verbose_lvl=verbose_lvl,
    trees = trees,
    dim = ncol(data),
    weights = weights,
    indexes = indexes,
    pmf = pmf,
    norm_matrix = norm_matrix,
    oob_pmf = oob_pmf,
    oob_kl = oob_kl,
    oob_ise = oob_ise
  ))
}

setMethod(f = "show", signature = c(object = "CortForest"), definition = function(object){
  cat(paste0("CortForest copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
             length(object@trees)," Cort trees."))
})

#' @describeIn rCopula-methods Method for the class CortForest
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "CortForest"), definition = function(n, copula) {
  if(n==0){
    return(matrix(0,nrow=0,ncol=copula@dim))
  }
  sampled_indexes = resample(1:length(copula@trees),size=n,prob = copula@weights,replace=TRUE)
  return(do.call(rbind,purrr::map(unique(sampled_indexes),~rCopula(sum(sampled_indexes==.x),copula@trees[[.x]]))))

})

#' @describeIn pCopula-methods Method for the class CortForest
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "CortForest"), definition = function(u, copula) {
  return(as.vector(vapply(copula@trees,function(t){pCopula(u,t)},rep(0.5,nrow(u))) %*% copula@weights))
})

#' @describeIn dCopula-methods Method for the class CortForest
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="CortForest"),   definition = function(u, copula) {
  return(vapply(copula@trees,function(t){dCopula(u,t)},rep(0.5,nrow(u))) %*% copula@weights)
})

#' @describeIn constraint_infl-methods Method for the class Cort
setMethod(f = "constraint_infl", signature = c(object="CortForest"),   definition = function(object) {
  return(purrr::map_dbl(object@trees,constraint_infl))
})

#' @describeIn quad_norm-methods Method for the class Cort
setMethod(f = "quad_norm", signature = c(object="CortForest"),   definition = function(object) {
  return(object@weights %*% object@norm_matrix %*% object@weights)
})









