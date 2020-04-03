#' @include generics.R empiricalCopula.R
NULL

.Cort = setClass(Class = "Cort", contains = "empiricalCopula",
                        slots = c(p_value_for_dim_red = "numeric",
                                  number_max_dim = "numeric",
                                  min_node_size = "numeric",
                                  verbose_lvl="numeric",
                                  #leaves = "list",
                                  vols = "numeric",
                                  f = "numeric",
                                  p = "numeric",
                                  a = "matrix",
                                  b = "matrix"),
                        validity = function(object) {
                          errors <- c()
                          if(!(length(dim(object@data))==2)){
                            errors <- c(errors, "data should be a matrix")
                          }
                          if(!(all(object@data >= 0) & all(1 > object@data))){
                            errors <- c(errors,"data should be a matrix inside 0,1")
                          }
                          if(object@number_max_dim > ncol(object@data)){
                            errors <- c(errors, "the maximum number of splitting dimensions should be smaller than the umber of dimensons!")
                          }
                          if(object@number_max_dim < 2){
                            errors <- c(errros, "splits cannot be done in dimmensions smaller than 2.")
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
#' @param slsqp_options options for nloptr::slsqp to find breakpoints : you can change defaults.
#' @param verbose_lvl numeric. set the verbosity. 0 for no ouptut and bigger you set it the most output you get.
#'
#' @name Cort-Class
#' @title The Cort estimator
#' @rdname Cort-Class
#'
#' @return a Cort object that can be fitted easily to produce a copula estimate.
#' @export
#'
#' @examples
#' (Cort(LifeCycleSavings[,1:3]))
Cort = function(x,
                p_value_for_dim_red=0.75,
                min_node_size=1,
                pseudo_data=FALSE,
                number_max_dim=NULL,
                verbose_lvl=1,
                slsqp_options = NULL) {

  # coerce the data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL

  if(!pseudo_data){
    data = apply(data,2,rank,ties.method="random")/(nrow(data)+1)
  }

  d = ncol(data)
  #construct the main leave :
  model = .Cort(
    data = data,
    p_value_for_dim_red = p_value_for_dim_red,
    number_max_dim = min(number_max_dim,d),
    min_node_size = min_node_size,
    verbose_lvl=verbose_lvl,
    dim = ncol(data),
    vols = c(1),
    f = c(1),
    p = c(1),
    a = matrix(0,ncol=d,nrow=1),
    b = matrix(1,ncol=d,nrow=1)
  )

  # then fit and return it :
  return(fit(model,slsqp_options))
}

setMethod(f = "show", signature = c(object = "Cort"), definition = function(object){
            cat(paste0("Cort copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
                       nrow(object@a)," leaves."))
          })

setMethod(f="fit", signature = c(object="Cort"), definition = function(object,slsqp_options = NULL, N = 999){


            # Deal with slsqp options :
            DEFAULT_SLQP_OPTIONS = list(stopval = -Inf,xtol_rel = 1e-4,maxeval = 100000,ftol_rel = 1e-6,ftol_abs = 1e-6)
            # get setted options :
            if(is.null(slsqp_options)){ slsqp_options = DEFAULT_SLQP_OPTIONS}
            if(!is.null(slsqp_options)){
              if(is.null(slsqp_options$stopval))  slsqp_options$stopval  = DEFAULT_SLQP_OPTIONS$stopval
              if(is.null(slsqp_options$xtol_rel)) slsqp_options$xtol_rel = DEFAULT_SLQP_OPTIONS$xtol_rel
              if(is.null(slsqp_options$maxeval))  slsqp_options$maxeval  = DEFAULT_SLQP_OPTIONS$maxeval
              if(is.null(slsqp_options$ftol_rel)) slsqp_options$ftol_rel = DEFAULT_SLQP_OPTIONS$ftol_rel
              if(is.null(slsqp_options$ftol_abs)) slsqp_options$ftol_abs = DEFAULT_SLQP_OPTIONS$ftol_abs
            }



            # Splitting the domain into small boxes :
            if(object@verbose_lvl>0) {cat("Splitting...\n")}

            # initialise the first leave [a,b]
            d = object@dim
            object@a = matrix(0,nrow=1,ncol=d)
            object@b = matrix(1,nrow=1,ncol=d)
            data_dist = list(object@data)
            split_dims = list(1:d)

            # splitting loop :
            continue = TRUE
            while(continue){

              are_splittables = purrr::map_lgl(1:nrow(object@a),function(j){
                (nrow(data_dist[[j]])>object@min_node_size) && (length(split_dims[[j]])>1)
              })

              if(!any(are_splittables)){
                continue = FALSE # we get out the while loop
              } else {
                if(object@verbose_lvl>0){cat("\n    ",sum(are_splittables),"leaves to split...")}
                if(object@verbose_lvl>1){cat("\n")}
                i_leaf_to_remove = numeric()
                for(i_leaf in which(are_splittables)){

                  # verbosity :
                  if(object@verbose_lvl>1){
                    cat(paste0("        Leaf with ",nrow(data_dist[[i_leaf]])," points.\n"))
                  }
                  if(object@verbose_lvl>2){
                    verb_df = data.frame(min = object@a[i_leaf,], max = object@b[i_leaf,])
                    verb_df$bp      = rep(NaN,d)
                    verb_df$p_value = rep(NaN,d)
                    verb_df$action  = rep("",d)
                    verb_df$reason  = rep("",d)
                    row.names(verb_df) = paste0("             ",1:nrow(verb_df))
                  }

                  # Randomize splitting dimensions :
                  if(length(split_dims[[i_leaf]]) > object@number_max_dim){
                    random_dims       = sample(x =split_dims[[i_leaf]], size=object@number_max_dim, replace=FALSE)
                    non_taken_dims    = split_dims[[i_leaf]][!(split_dims[[i_leaf]] %in% random_dims)]
                    split_dims[[i_leaf]] = random_dims
                  } else{
                    non_taken_dims = numeric()
                  }

                  # verbosity :
                  if(object@verbose_lvl>2){
                    verb_df$action[non_taken_dims] = "Dissmissed"
                    verb_df$reason[non_taken_dims] = "Randomly"
                  }

                  if(length(split_dims[[i_leaf]])<=1){
                    # we just keep the leaf :
                    split_dims[[i_leaf]] = c(split_dims[[i_leaf]],non_taken_dims)
                  } else { # we try to split:

                    # Compute prerequisites for the optimisation of the breakpoint
                    n = nrow(data_dist[[i_leaf]])
                    d_split = length(split_dims[[i_leaf]])
                    a = object@a[i_leaf,split_dims[[i_leaf]]]
                    b = object@b[i_leaf,split_dims[[i_leaf]]]
                    z = (t(data_dist[[i_leaf]][,split_dims[[i_leaf]]]) - a)/(b-a) # d*n
                    bin_repr = sapply(1:(2^d_split),number2binary,d_split)

                    #Launche optimisation routine for the breakpoint :
                    optimizer = nloptr::slsqp(
                      x0 = rowMeans(z),
                      fn = lossFunc, # a cpp function
                      lower = rep(0,d_split),
                      upper = rep(1,d_split),
                      nl.info=object@verbose_lvl>3,
                      control=slsqp_options,
                      bin_repr = bin_repr,
                      z = z)

                    # Get the breakpoints and the final splitting :
                    bp = a + optimizer$par*(b-a)
                    z_min = bp*bin_repr
                    z_max = bp^(1-bin_repr)

                    # Compute p-values for the breakpoint :
                    montecarlo = cortMonteCarlo(z,z_min,z_max,as.integer(N)) # a cpp function
                    p_values = rowMeans(montecarlo[1,] <= t(montecarlo[-1,]))

                    if(any(is.na(p_values))){
                      p_values[is.na(p_values)] = 0
                    }

                    if(object@verbose_lvl>2){
                      verb_df$bp[split_dims[[i_leaf]]]      <- bp
                      verb_df$p_value[split_dims[[i_leaf]]] <- p_values
                    }

                    # if p_values are too big, remove dimensions.
                    p_val_too_big  = p_values > object@p_value_for_dim_red
                    if(object@verbose_lvl>2){
                      verb_df$action[split_dims[[i_leaf]][p_val_too_big]] = "Removed"
                      verb_df$reason[split_dims[[i_leaf]][p_val_too_big]] = "Independence test"
                    }

                    # if any of the breakpoints are too close to boundary, we remove the dimensions :
                    normed_bp      = (bp - object@a[i_leaf,split_dims[[i_leaf]]])/(object@b[i_leaf,split_dims[[i_leaf]]] - object@a[i_leaf,split_dims[[i_leaf]]])
                    threshold          = 1/min((nrow(data_dist[[i_leaf]])+1)^2,1000)
                    close_to_bound = (normed_bp< threshold) + (normed_bp > 1-threshold) > 0

                    if(object@verbose_lvl>2){
                      verb_df$action[split_dims[[i_leaf]][close_to_bound]] = "Removed"
                      verb_df$reason[split_dims[[i_leaf]][close_to_bound]] = "Close to boundary"
                    }

                    to_be_removed = p_val_too_big+close_to_bound>0

                    if(all(to_be_removed)){
                      # we just keep the leaf :
                      split_dims[[i_leaf]] = non_taken_dims
                    } else {

                      if(any(to_be_removed)){
                        split_dims[[i_leaf]] = split_dims[[i_leaf]][!to_be_removed]
                        bp = bp[!to_be_removed]
                      }

                      if(length(split_dims[[i_leaf]]) == 1){
                        # we just keep the leaf as is :
                        if(object@verbose_lvl>2){
                          verb_df$action[split_dims[[i_leaf]]] = "Dissmissed"
                          verb_df$reason[split_dims[[i_leaf]]] = "No one-dim split"
                        }
                        split_dims[[i_leaf]] = c(split_dims[[i_leaf]],non_taken_dims)
                      } else {
                        # NOW WE SPLIT
                        i_leaf_to_remove = c(i_leaf_to_remove,i_leaf)
                        if(object@verbose_lvl>2){verb_df$action[split_dims[[i_leaf]]] = "Splitted"}

                        # remove the breakpoint from the data points if it's one of them :
                        are_the_breakpoint  = (colSums(t(data_dist[[i_leaf]][,split_dims[[i_leaf]]]) == bp) == length(split_dims[[i_leaf]]))
                        if(any(are_the_breakpoint)){
                          if(object@verbose_lvl>4){cat("be carrefull, we are splitting on a point.\n")}
                          data_dist[[i_leaf]] = data_dist[[i_leaf]][!are_the_breakpoint,]
                        }

                        # construct new information for new leaves :
                        d_split=length(split_dims[[i_leaf]])
                        D = 2^d_split
                        for (i in 1:D){
                          bin = number2binary(i,d_split)
                          new_a = object@a[i_leaf,]
                          new_b = object@b[i_leaf,]
                          new_a[split_dims[[i_leaf]]] = bp^bin * new_a[split_dims[[i_leaf]]] ^ (1-bin)
                          new_b[split_dims[[i_leaf]]] = bp^(1-bin) * new_b[split_dims[[i_leaf]]] ^ bin
                          i_new_data = apply((t(data_dist[[i_leaf]]) >= new_a)*(t(data_dist[[i_leaf]])<new_b),2,all)
                          new_data = data_dist[[i_leaf]][i_new_data, , drop=FALSE]

                          # imput the new information from the leaf :
                          object@a = rbind(object@a,new_a)
                          object@b = rbind(object@b,new_b)
                          split_dims = c(split_dims,list(c(split_dims[[i_leaf]],non_taken_dims))) # append new split dims
                          data_dist = c(data_dist,list(new_data))
                        }
                      }
                    }
                  }

                  if(object@verbose_lvl>2) {
                    cat(toString.data.frame(verb_df,digits=8))
                    cat("\n")
                  }
                  if(object@verbose_lvl>2){cat("\n")}

                }
                # remove information from the splitted leaves :
                if(length(i_leaf_to_remove) >= 1){
                  object@a = object@a[-i_leaf_to_remove,]
                  object@b = object@b[-i_leaf_to_remove,]
                  split_dims = split_dims[-i_leaf_to_remove]
                  data_dist = data_dist[-i_leaf_to_remove]
                }
              }
            }
            #  info about the data distribution :
            object@vols = apply(object@b-object@a,1,prod)
            object@f    = purrr::map_dbl(data_dist,~nrow(.x))/nrow(object@data)

            # Then compute the weights :
            if(object@verbose_lvl>0) {cat("\nEnforcing constraints...\n")}

            # Building constraints :
            n = nrow(object@a)
            d = ncol(object@a)
            evaluation_points = purrr::map(1:d,~unique((object@b+object@a)[,.x]/2))
            dims = unlist(purrr::map2(evaluation_points,1:d,function(x,y){rep(y,length(x))}))
            F_vec = unlist(evaluation_points)
            lambdas = pmin(pmax((F_vec - t(object@a[,dims,drop=FALSE]))/(t(object@b[,dims,drop=FALSE] - object@a[,dims,drop=FALSE])),0),1)

            # Constructing the parameters for the osqp solver :
            P_mat = diag(1/object@vols)
            q_vec = -object@f/object@vols
            A_mat = rbind(lambdas,rep(1,n),diag(n))
            l_vec = c(F_vec,1,rep(0,n))
            u_vec = c(F_vec,1,rep(Inf,n))

            # building the model
            if (object@verbose_lvl>1) { model = osqp::osqp(P=P_mat, q=q_vec, A=A_mat, l=l_vec, u=u_vec, pars=osqp::osqpSettings(max_iter = 100000L,
                                 eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = TRUE))
            } else { model = osqp::osqp(P=P_mat, q=q_vec, A=A_mat, l=l_vec, u=u_vec, pars=osqp::osqpSettings(max_iter = 100000L,
                                 eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = FALSE))
            }

            # Launching the solver
            model$WarmStart(x=object@f)
            rez = model$Solve()
            # saving weights:
            object@p = pmax(rez$x,0)/sum(pmax(rez$x,0)) # correction for small negative weights.

            if(object@verbose_lvl>0){cat("Done !\n")}
            return(object)
          })

#' @describeIn rCopula-methods Method for the class Cort
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "Cort"), definition = function(n, copula) {

  # We need to simulate n rnadom vctors from the fitted model :
  # for that, we will samples indexes of the boxes according to weights :
  sampled_indexes = sample(1:nrow(copula@a),size=n,prob = copula@p,replace=TRUE)
  sim = matrix(runif(n*copula@dim),nrow=n,ncol=copula@dim)
  sim = copula@a[sampled_indexes,] + sim * (copula@b[sampled_indexes,]-copula@a[sampled_indexes,])
  return(sim)
})

#' @describeIn pCopula-methods Method for the class Cort
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "Cort"), definition = function(u, copula) {

  # for the c.d.f of the copula, we need to compute the measure_in per box and sum them :
  a = copula@a
  b = copula@b
  p = copula@p

  n = dim(a)[1]
  d = dim(a)[2]
  m = dim(u)[1]

  dim(u) = c(m,1,d)
  dim(a) = c(1,n,d)
  dim(b) = c(1,n,d)
  dim(p) = c(1,n)

  core = (u[,rep(1,n),,drop=FALSE] - a[rep(1,m),,,drop=FALSE])/(b-a)[rep(1,m),,,drop=FALSE] # m, n, d
  measure_in = apply(pmax(pmin(core,1),0),1:2,prod) # m,n
  return(rowSums(measure_in * p[rep(1,m),,drop=FALSE])) # m

  #return(rowSums(vapply(copula@leaves,function(l){measure_in(l,u)*l@weight},FUN.VALUE=numeric(nrow(u)))))
})

#' @describeIn dCopula-methods Method for the class Cort
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="Cort"),   definition = function(u, copula) {

  a = copula@a
  b = copula@b
  p = copula@p
  vols = copula@vols

  n = dim(a)[1]
  d = dim(a)[2]
  m = dim(u)[1]

  dim(u) = c(m,1,d)
  dim(a) = c(1,n,d)
  dim(b) = c(1,n,d)
  dim(p) = c(1,n)
  dim(vols) = c(1,n)

  core = (u[,rep(1,n),,drop=FALSE] >= a[rep(1,m),,,drop=FALSE])*(u[,rep(1,n),,drop=FALSE] < b[rep(1,m),,,drop=FALSE]) # m, n, d
  are_in = (rowSums(core,dims = 2)==d)
  return(rowSums(are_in * p[rep(1,m),,drop=FALSE]/vols[rep(1,m),,drop=FALSE])) # m




 # return(rowSums(vapply(copula@leaves,function(l){contains(l,u,type="rllc")*l@weight},FUN.VALUE = rep(0.5,nrow(u)))))

})

#' @describeIn biv_rho-methods Method for the class Cort
setMethod(f = "biv_rho", signature = c(copula="Cort"),   definition = function(copula) {
  # should return the spearmann rho for this copula.

  d = ncol(copula@data)
  n = nrow(copula@a)
  rho = diag(d)
  p = copula@p
  a = copula@a
  b = copula@b

  for (i in 2:d){
    for (j in 1:(i-1)){
      rho[j,i] <- rho[i,j] <- 12*sum(apply((2-b[,c(i,j),drop=FALSE] - a[,c(i,j),drop=FALSE])/2,1,prod)*p)-3
    }
  }

  return(rho)


})

# helper function for the kendall tau computation :
B <- function(a,b,n_leaves,d_dim){

  c = a
  d = b

  dim(a) <- c(n_leaves, 1,        d_dim)
  dim(b) <- c(n_leaves, 1,        d_dim)
  dim(c) <- c(1,        n_leaves, d_dim)
  dim(d) <- c(1,        n_leaves, d_dim)

  a <- a[,                rep(1,n_leaves),,drop=FALSE]
  b <- b[,                rep(1,n_leaves),,drop=FALSE]
  c <- c[rep(1,n_leaves),                ,,drop=FALSE]
  d <- d[rep(1,n_leaves),                ,,drop=FALSE]

  x = pmax(a,c)
  y = pmin(b,d)
  z = pmax(a,d)

  return((y *( y/2 - c)  - x * (x/2 -c)) * (x < y) + (b - z) * (d - c) * (z < b))
}

#' @describeIn biv_tau-methods Method for the class Cort
setMethod(f = "biv_tau", signature = c(copula="Cort"),   definition = function(copula) {
  # dimensions :
  d = ncol(copula@data)
  n = nrow(copula@a)
  dims = cbind(rep(1:d,d),rep(1:d,each=d))
  dims = dims[dims[,1]>dims[,2],,drop=FALSE]
  K = nrow(dims)

  # Extract info from the model :
  weights = copula@p
  a = copula@a
  b = copula@b

  # comput ethe cross_kernel :
  cross_B = B(a,b,n,d) # n, n, d
  cross_B = apply(vapply(1:K,function(i){cross_B[,,dims[i,],drop=FALSE]},cross_B[,,dims[1,],drop=FALSE]),c(1,2,4),prod) # n,n,K

  # compute boxes measures :
  side_lengths = b - a #n, d
  measures = t(apply(vapply(1:K,function(i){side_lengths[,dims[i,],drop=FALSE]},side_lengths[,c(1,2),drop=FALSE]),c(1,3),prod)) # K, n

  # compute the cross_kernel of the model :
  kernel = weights/t(measures) # n,K
  kernel2 = kernel
  dim(kernel) = c(n,1,K)
  dim(kernel2) = c(1,n,K)

  # finaly :
  cross_kernel = aperm(kernel[,rep(1,n),,drop=FALSE]*kernel2[rep(1,n),,,drop=FALSE]*cross_B,c(3,1,2)) # K, n, n
  kappa = 4 * rowSums(cross_kernel) - 1

  # just put them back in a 2*2 matrix :
  mat = diag(d)/2
  mat[dims] <- kappa
  mat = mat+t(mat)
  return(mat)

})

#' @describeIn loss-methods Method for the class Cort
setMethod(f = "loss", signature = c(object="Cort"),   definition = function(object) {
  return(sum((object@p^2/2 - object@p*object@f)/object@vols))
})

#' @describeIn constraint_infl-methods Method for the class Cort
setMethod(f = "constraint_infl", signature = c(object="Cort"),   definition = function(object) {
  return(sum(((object@p - object@f)^2/object@vols)))
})

#' @describeIn quad_norm-methods Method for the class Cort
setMethod(f = "quad_norm", signature = c(object="Cort"),   definition = function(object) {
  return(sum(object@p^2/object@vols))
})

#' @describeIn quad_prod_with_data-methods Method for the class Cort
setMethod(f = "quad_prod_with_data", signature = c(object="Cort"),   definition = function(object) {
  return(sum(object@f * object@p/object@vols))
})

#' @describeIn quad_prod-methods Method for the class Cort
setMethod(f = "quad_prod", signature = c(object="Cort",other_tree = "Cort"),   definition = function(object,other_tree) {

  # The implemntation is in C++
  return(quadProd(
    object@a,
    object@b,
    object@p/object@vols,
    other_tree@a,
    other_tree@b,
    other_tree@p/other_tree@vols
  ))

})

#' @param M the number of simulations
#'
#' @describeIn kendall_func-methods Method for the class Cort
setMethod(f = "kendall_func", signature = c(object="Cort"),   definition = function(object,t,M=1000) {

  # Simulate a bunch of random variables :
  m = length(t)
  rng = pCopula(rCopula(M,object),object)
  dim(rng) = c(M,1)
  dim(t) = c(1,m)
  colMeans(rng[,rep(1,m),drop=FALSE]<=t[rep(1,M),,drop=FALSE])

})

#' @describeIn project_on_dims-methods Method for the class Cort
setMethod(f = "project_on_dims", signature = c(object="Cort"),   definition = function(object,dims) {

  # this function should project the tree on a smaller subset of dimensions.
  # for the moment only sets of 2 dimensions are supported, although bigger sets might be used;
  if(length(dims) != 2){
    stop("only two-dimensional projection are supported for the moment")
  }

  # first, getting the edges of the bins :
  a       = object@a
  b       = object@b
  p_over_vol = object@p/object@vols
  n = nrow(a)
  d = ncol(a)

  points = rbind(a,b)
  edges = purrr::map(dims,function(i){
    sort(unique(points[,i]))
  }) # a list; dimension by dimensions, of the breakpoints in the dimensions we want to keep.

  # now we need ot construct the boxes :
  ks = 1:(length(edges[[1]])-1)
  ls = 1:(length(edges[[2]])-1)

  new_n = length(ks)*length(ls)
  new_a = matrix(0,nrow=new_n,ncol=2)
  new_b = new_a
  i=1
  for (k in ks){
    for (l in ls){
      new_a[i,] <- c(edges[[1]][k  ], edges[[2]][l  ])
      new_b[i,] <- c(edges[[1]][k+1], edges[[2]][l+1])
      i = i+1
    }
  }

  a_complete = rep(0,d)
  b_complete = rep(1,d)

  # set things :
  object@data = object@data[,dims]
  object@f = purrr::map_dbl(1:new_n,function(i){
    sum(colSums((new_a[i,] <= t(object@data)) & (t(object@data) < new_b[i,]))==2)
  })/nrow(object@data)
  object@p = purrr::map_dbl(1:new_n,function(i){
      a_complete[dims] <- new_a[i,]
      b_complete[dims] <- new_b[i,]
      return(sum(apply(pmax(pmin(t(b),b_complete)-pmax(t(a),a_complete),0),2,prod)*p_over_vol))
    })
  object@dim = 2
  object@a = new_a
  object@b = new_b
  object@vols = apply((new_b - new_a),1,prod)

  return(object)
})

#' @export
plot.Cort <- function(x,...){

  d = ncol(x@data)
  dd = do.call(rbind,unlist(purrr::map(1:(d-1),function(i){
    purrr::map((i+1):(d),function(j){
      proj = project_on_dims(x,c(i,j))
      rbind(
        data.frame(xmin = proj@a[,1],
                   ymin = proj@a[,2],
                   xmax = proj@b[,1],
                   ymax = proj@b[,2],
                   weight = proj@p,
                   volume = proj@vols,
                   dim_x = rep(i,nrow(proj@a)),
                   dim_y = rep(j,nrow(proj@a))),
        data.frame(ymin = proj@a[,1],
                   xmin = proj@a[,2],
                   ymax = proj@b[,1],
                   xmax = proj@b[,2],
                   weight = proj@p,
                   volume = proj@vols,
                   dim_y = rep(i,nrow(proj@a)),
                   dim_x = rep(j,nrow(proj@a)))
      )
    })}),recursive=FALSE))
    dd$col = log(1+dd$weight/dd$volume)
    dd$col = dd$col/max(dd$col)

    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    graphics::par(mfrow=c(d,d),
      mai = c(0,0,0,0),
      oma = c(3,3,5,3))

  for (i in 1:d){
    for (j in 1:d){
      graphics::plot(c(0,1),c(0,1),type="n",xlab="",ylab="",xaxt='n',yaxt='n')
      if(i != j){
        xx = dd[(dd$dim_x == i)&(dd$dim_y==j),]
        graphics::rect(xx$ymin,xx$xmin,xx$ymax,xx$xmax,col=grDevices::gray(1-xx$col),border=NA,density=NA)

          graphics::points(x@data[,j],x@data[,i],cex=0.5,col="red")

      }
    }
  }
}





























