#' @include generics.R Box.R WeightedBox.R
NULL

.Cort = setClass(Class = "Cort",
                        slots = c(data = "matrix",
                                  p_value_for_dim_red = "numeric",
                                  number_max_dim = "numeric",
                                  verbose="logical",
                                  leaves = "list"),
                        validity = function(object) {
                          errors <- c()
                          if(!(length(dim(object@data))==2)){
                            errors <- c(errors, "data should be a matrix")
                          }
                          if(!(all(object@data >= 0) & all(1 > object@data))){
                            errors <- c(errors,"data should be a matrix inside 0,1")
                          }
                          if(object@number_max_dim > ncol(object@data)){
                            errors <- c(errors, "the maimum number of splitting dimensions should be smaller than the umber of dimensons!")
                          }
                          if (length(errors) == 0)
                            TRUE else errors
                        })

#' Cort class
#'
#' This class implements the CORT algorithm to a fit a multivariate copula using piece constant density.
#'
#'
#' @param x The data, must be provided as a matrix with each row as an observation.
#' @param p_value_for_dim_red a p_value for the localised dimension reduction test
#' @param min_node_size The minimum number of observation avaliable in a leaf to initialise a split.
#' @param pseudo_data set to True if you are already providing data on the copula space.
#' @param number_max_dim The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.
#' @param verbose set to TRUE for more detailed output during the fit.
#' @name Cort-Class
#' @title The Cort estimator
#' @rdname Cort-Class
#'
#' @return a Cort object that can be fitted easily to produce a copula estimate.
#' @export
#'
#' @examples
#' library(cort)
#' (Cort(LifeCycleSavings[,1:3]))
Cort = function(x,
                p_value_for_dim_red=0.75,
                min_node_size=1,
                pseudo_data=FALSE,
                number_max_dim=NULL,
                verbose=TRUE) {

  # coerce the data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL

  if(!pseudo_data){
    data = apply(data,2,rank)/(nrow(data)+1)
  }

  d = ncol(data)
  #construct the main leave :
  root = WeightedBox(Box(rep(0,d),rep(1,d)),data)
  model = .Cort(
    data = data,
    p_value_for_dim_red = p_value_for_dim_red,
    number_max_dim = max(number_max_dim,d),
    leaves=list(root),
    verbose=verbose
  )

  # then fit and return it :
  return(fit(model))
}

#' @describeIn Cort-Class dimension of the copula
setMethod(f = "dim", signature = (x = "Cort"), definition = function(x){return(ncol(x@data))})

setMethod(f = "show",
          signature = c(object = "Cort"),
          definition = function(object){
            cat(paste0("Cort copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
                       length(object@leaves)," leaves"))
          })

setMethod(f="fit",
          signature = c(object="Cort"),
          definition = function(object){

            # first, split the boxes:
            if(object@verbose) {cat("Splitting...\n")}
            continue = TRUE
            while(continue){
              are_splittables = purrr::map_lgl(object@leaves,is_splittable)
              if(object@verbose){cat("    ",sum(are_splittables),"leaves to split...\n")}
              if(any(are_splittables)){

                # Split every box :
                object@leaves = purrr::map(object@leaves,
                                           ~split(.x,
                                                  p_val_threshold = object@p_value_for_dim_red,
                                                  number_max_dim=object@number_max_dim,
                                                  verbose=object@verbose
                  )) %>%
                  unlist(recursive = FALSE)


                continue = TRUE
              } else {
                continue = FALSE
              }
            }

            # Then compute the weights :
            if(object@verbose) {cat("Enforcing constraints...\n")}

            #  Objective :

            volumes = purrr::map_dbl(object@leaves,"volume")
            starting_point = purrr::map_dbl(object@leaves,"weight")
            contain_info = sapply(object@leaves,function(l){contains(l,object@data)})
            P_mat = diag(1/volumes)
            q_vec = -rowMeans(t(contain_info)/volumes)

            # Constraints :
            n = length(object@leaves)
            as = sapply(object@leaves,function(x){x@a})
            bs = sapply(object@leaves,function(x){x@b})
            middle_points = (bs+as)/2
            evaluation_points = apply(middle_points,1,unique)

            dims = purrr::map2(evaluation_points,1:ncol(object@data),function(x,y){rep(y,length(x))}) %>% unlist()

            F_vec = unlist(evaluation_points)

            lambdas = pmin(pmax((F_vec - as[dims,])/(bs[dims,] - as[dims,]),0),1)

            A_mat = rbind(lambdas,rep(1,n),diag(n))

            l_vec = c(F_vec,1,rep(0,n))
            u_vec = c(F_vec,1,rep(Inf,n))

            model = build_model(P_mat,q_vec,A_mat,l_vec,u_vec,object@verbose)
            model$WarmStart(x=starting_point)
            rez = model$Solve()

            new_weights = rez$x
            new_weights = pmax(new_weights,0)/sum(pmax(new_weights,0)) # correction for small negative weights.

            object@leaves = purrr::map2(object@leaves, new_weights,function(l,w){
              l@weight = w
              return(l)
            })


            if(object@verbose){cat("Done !\n")}
            return(object)
          })



build_model <- function(P_mat,q_vec,A_mat,l_vec,u_vec,verbose=TRUE){

  if(verbose){
    return(osqp::osqp(P=P_mat,
                      q=q_vec,
                      A=A_mat,
                      l=l_vec,
                      u=u_vec,
                      pars=osqp::osqpSettings(rho = 0.1,
                                              sigma = 1e-06,
                                              max_iter = 100000L,
                                              eps_abs = 0.000001,
                                              eps_rel = 0.000001,
                                              eps_prim_inf = 0.000001,
                                              eps_dual_inf = 0.000001,
                                              alpha = 1.6,
                                              linsys_solver = 0L, #0 for qdldl, 1 for mkl_paradiso
                                              delta = 1e-06,
                                              polish = TRUE,
                                              polish_refine_iter = 100L,
                                              verbose = TRUE,
                                              scaled_termination = FALSE,
                                              check_termination = 25L,
                                              warm_start = TRUE,
                                              scaling = 10L,
                                              adaptive_rho = 1L,
                                              adaptive_rho_interval = 0L,
                                              adaptive_rho_tolerance = 5,
                                              adaptive_rho_fraction = 0.4)))
  }
  if(!verbose){
    return(osqp::osqp(P=P_mat,
                      q=q_vec,
                      A=A_mat,
                      l=l_vec,
                      u=u_vec,
                      pars=osqp::osqpSettings(rho = 0.1,
                                              sigma = 1e-06,
                                              max_iter = 100000L,
                                              eps_abs = 0.000001,
                                              eps_rel = 0.000001,
                                              eps_prim_inf = 0.000001,
                                              eps_dual_inf = 0.000001,
                                              alpha = 1.6,
                                              linsys_solver = 0L, #0 for qdldl, 1 for mkl_paradiso
                                              delta = 1e-06,
                                              polish = TRUE,
                                              polish_refine_iter = 100L,
                                              verbose = FALSE,
                                              scaled_termination = FALSE,
                                              check_termination = 25L,
                                              warm_start = TRUE,
                                              scaling = 10L,
                                              adaptive_rho = 1L,
                                              adaptive_rho_interval = 0L,
                                              adaptive_rho_tolerance = 5,
                                              adaptive_rho_fraction = 0.4)))
  }
}








































