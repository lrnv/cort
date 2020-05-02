#' @include generics.R empiricalCopula.R
NULL

.Cort = setClass(Class = "Cort", contains = "empiricalCopula",
                        slots = c(p_value_for_dim_red = "numeric",
                                  number_max_dim = "numeric",
                                  min_node_size = "numeric",
                                  verbose_lvl="numeric",
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
                            errors <- c(errors, "the maximum number of splitting dimensions should be smaller than the number of dimensons!")
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
#' @param N The number of bootstrap resamples for p_values computations.
#' @param osqp_options options for the weights optimisation. You can pass a call to osqp::osqpSettings, or NULL for defaults.
#' @param force_grid boolean. set to TRUE to force breakpoint to be on the n-checkerboard grid.
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
                slsqp_options = NULL,
                osqp_options = NULL,
                N = 999,
                force_grid=FALSE) {

  # Coerce the data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL

  if(!pseudo_data){
    data = apply(data,2,rank,ties.method="random")/(nrow(data)+1)
  }

  d = ncol(data)
  n_obs = nrow(data)

  # Construct the object :
  object = .Cort(
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

  # Deal with solver parameters :
  DEFAULT_SLQP_OPTIONS = list(stopval = -Inf,xtol_rel = 1e-4,maxeval = 100000,ftol_rel = 1e-6,ftol_abs = 1e-6)
  if(is.null(slsqp_options)){ slsqp_options = DEFAULT_SLQP_OPTIONS}
  if(!is.null(slsqp_options)){
    if(is.null(slsqp_options$stopval))  slsqp_options$stopval  = DEFAULT_SLQP_OPTIONS$stopval
    if(is.null(slsqp_options$xtol_rel)) slsqp_options$xtol_rel = DEFAULT_SLQP_OPTIONS$xtol_rel
    if(is.null(slsqp_options$maxeval))  slsqp_options$maxeval  = DEFAULT_SLQP_OPTIONS$maxeval
    if(is.null(slsqp_options$ftol_rel)) slsqp_options$ftol_rel = DEFAULT_SLQP_OPTIONS$ftol_rel
    if(is.null(slsqp_options$ftol_abs)) slsqp_options$ftol_abs = DEFAULT_SLQP_OPTIONS$ftol_abs
  }
  if(is.null(osqp_options)){
    if (object@verbose_lvl>1) {
      osqp_options = osqp::osqpSettings(max_iter = 100000L, eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = TRUE)
    } else {
      osqp_options = osqp::osqpSettings(max_iter = 100000L, eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = FALSE)
    }
  }

  # Now we start fitting.
  if(object@verbose_lvl>0) {cat("Splitting...\n")}
  dd = list(object@data)
  splitd = list(1:d) # the splitting dimensions

  # Loop until there are no more splittable leaves :
  continue = TRUE
  while(continue){

    are_splittables = purrr::map_lgl(1:nrow(object@a),function(j){
      (nrow(dd[[j]])>object@min_node_size) && (length(splitd[[j]])>1)
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
          cat(paste0("        Leaf with ",nrow(dd[[i_leaf]])," points.\n"))
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
        if(length(splitd[[i_leaf]]) > object@number_max_dim){
          random_dims       = resample(x =splitd[[i_leaf]], size=object@number_max_dim, replace=FALSE)
          non_taken_dims    = splitd[[i_leaf]][!(splitd[[i_leaf]] %in% random_dims)]
          splitd[[i_leaf]] = random_dims
        } else{
          non_taken_dims = numeric()
        }

        # verbosity :
        if(object@verbose_lvl>2){
          verb_df$action[non_taken_dims] = "Dissmissed"
          verb_df$reason[non_taken_dims] = "Randomly"
        }

        if(length(splitd[[i_leaf]])<=1){
          # we just keep the leaf :
          splitd[[i_leaf]] = c(splitd[[i_leaf]],non_taken_dims)
        } else { # we try to split:

          # Compute prerequisites for the optimisation of the breakpoint
          n = nrow(dd[[i_leaf]])
          d_split = length(splitd[[i_leaf]])
          a = object@a[i_leaf,splitd[[i_leaf]]]
          b = object@b[i_leaf,splitd[[i_leaf]]]
          z = (t(dd[[i_leaf]][,splitd[[i_leaf]]]) - a)/(b-a) # d*n
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
          z_bp = optimizer$par
          bp = a + z_bp*(b-a)

          if(force_grid){
            bp = round(bp*2*(n_obs+1))/(n_obs+1)/2
            z_bp = (bp-a)/(b-a)
          }

          # Compute p-values for the breakpoint :
          z_min = z_bp*bin_repr
          z_max = z_bp^(1-bin_repr)
          p_values = cortMonteCarlo(z,z_min,z_max,as.integer(N))

          if(any(is.na(p_values))){
            p_values[is.na(p_values)] = 0
          }

          if(object@verbose_lvl>2){
            verb_df$bp[splitd[[i_leaf]]]      <- bp
            verb_df$p_value[splitd[[i_leaf]]] <- p_values
          }

          # if p_values are too big, remove dimensions.
          p_val_too_big  = p_values > object@p_value_for_dim_red
          if(object@verbose_lvl>2){
            verb_df$action[splitd[[i_leaf]][p_val_too_big]] = "Removed"
            verb_df$reason[splitd[[i_leaf]][p_val_too_big]] = "Independence test"
          }

          # if any of the breakpoints are too close to boundary, we remove the dimensions :
          normed_bp      = (bp - object@a[i_leaf,splitd[[i_leaf]]])/(object@b[i_leaf,splitd[[i_leaf]]] - object@a[i_leaf,splitd[[i_leaf]]])
          threshold      = 1/ifelse(force_grid,(n_obs+1)^2,min((nrow(dd[[i_leaf]])+1)^2,1000))
          close_to_bound = (normed_bp< threshold) + (normed_bp > 1-threshold) > 0

          if(object@verbose_lvl>2){
            verb_df$action[splitd[[i_leaf]][close_to_bound]] = "Removed"
            verb_df$reason[splitd[[i_leaf]][close_to_bound]] = "Close to boundary"
          }

          to_be_removed = p_val_too_big+close_to_bound>0

          if(all(to_be_removed)){
            # we just keep the leaf :
            splitd[[i_leaf]] = non_taken_dims
          } else {

            if(any(to_be_removed)){
              splitd[[i_leaf]] = splitd[[i_leaf]][!to_be_removed]
              bp = bp[!to_be_removed]
            }

            if(length(splitd[[i_leaf]]) == 1){
              # we just keep the leaf as is :
              if(object@verbose_lvl>2){
                verb_df$action[splitd[[i_leaf]]] = "Dissmissed"
                verb_df$reason[splitd[[i_leaf]]] = "No one-dim split"
              }
              splitd[[i_leaf]] = c(splitd[[i_leaf]],non_taken_dims)
            } else {
              # NOW WE SPLIT
              i_leaf_to_remove = c(i_leaf_to_remove,i_leaf)
              if(object@verbose_lvl>2){verb_df$action[splitd[[i_leaf]]] = "Splitted"}

              # remove the breakpoint from the data points if it's one of them :
              are_the_breakpoint  = (colSums(t(dd[[i_leaf]][,splitd[[i_leaf]]]) == bp) == length(splitd[[i_leaf]]))
              if(any(are_the_breakpoint)){
                if(object@verbose_lvl>4){cat("be carrefull, we are splitting on a point.\n")}
                dd[[i_leaf]] = dd[[i_leaf]][!are_the_breakpoint,,drop=FALSE]
              }

              # construct new information for new leaves :
              d_split=length(splitd[[i_leaf]])
              D = 2^d_split
              for (i in 1:D){
                # build the leaf :
                bin = number2binary(i,d_split)
                new_a = object@a[i_leaf,]
                new_b = object@b[i_leaf,]
                new_a[splitd[[i_leaf]]] = bp^bin * new_a[splitd[[i_leaf]]] ^ (1-bin)
                new_b[splitd[[i_leaf]]] = bp^(1-bin) * new_b[splitd[[i_leaf]]] ^ bin
                i_new_data = colSums((t(dd[[i_leaf]]) >= new_a)*(t(dd[[i_leaf]])<new_b))==d
                new_data = dd[[i_leaf]][i_new_data, , drop=FALSE]

                # imput the new information from the leaf :
                object@a = rbind(object@a,new_a)
                object@b = rbind(object@b,new_b)
                splitd = c(splitd,list(c(splitd[[i_leaf]],non_taken_dims))) # append new split dims to non-taken ones.
                dd = c(dd,list(new_data))
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
        splitd = splitd[-i_leaf_to_remove]
        dd = dd[-i_leaf_to_remove]
      }
    }
  }
  #  info about the data distribution :
  object@vols = apply(object@b-object@a,1,prod)
  object@f    = purrr::map_dbl(dd,~nrow(.x))/nrow(object@data)

  # Then compute the weights :
  if(object@verbose_lvl>0) {cat("\nEnforcing constraints...\n")}

  # Building constraints for optimisation of the weights:
  n = nrow(object@a)
  d = ncol(object@a)
  evaluation_points = purrr::map(1:d,~unique((object@b+object@a)[,.x]/2))
  dims = unlist(purrr::map2(evaluation_points,1:d,function(x,y){rep(y,length(x))}))
  F_vec = unlist(evaluation_points)
  lambdas = pmin(pmax((F_vec - t(object@a[,dims,drop=FALSE]))/(t(object@b[,dims,drop=FALSE] - object@a[,dims,drop=FALSE])),0),1)

  # building the weights
  model = osqp::osqp(P=diag(1/object@vols),
                     q=-object@f/object@vols,
                     A=rbind(lambdas,rep(1,n),diag(n)),
                     l=c(F_vec,1,rep(0,n)),
                     u=c(F_vec,1,rep(Inf,n)),
                     pars=osqp_options)

  # Launching the solver
  model$WarmStart(x=object@f)
  rez = model$Solve()
  # saving weights:
  object@p = pmax(rez$x,0)/sum(pmax(rez$x,0)) # correction for small negative weights.

  # remove names :
  row.names(object@a) <- NULL
  row.names(object@b) <- NULL
  names(object@vols) <- NULL

  if(object@verbose_lvl>0){cat("Done !\n")}
  return(object)
}

setMethod(f = "show", signature = c(object = "Cort"), definition = function(object){
  cat(paste0("Cort copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
             nrow(object@a)," leaves."))
})


#' @describeIn rCopula-methods Method for the class Cort
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "Cort"), definition = function(n, copula) {
  # we sample the boxes and then sample from them.
  sampled_indexes = resample(1:nrow(copula@a),size=n,prob = copula@p,replace=TRUE)
  matrix(runif(n*copula@dim,min = as.vector(copula@a[sampled_indexes,]),max=as.vector(copula@b[sampled_indexes,])),nrow=n,ncol=copula@dim)
})

#' @describeIn pCopula-methods Method for the class Cort
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "Cort"), definition = function(u, copula) {
  # The implementation is in Rcpp.
  pCort(copula@a,
        copula@b,
        copula@p,
        u)
})

#' @describeIn dCopula-methods Method for the class Cort
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="Cort"),   definition = function(u, copula) {
  # The implementation is in Rcpp
  dCort(copula@a,
        copula@b,
        copula@p/copula@vols,
        u)
})

#' @describeIn biv_rho-methods Method for the class Cort
setMethod(f = "biv_rho", signature = c(copula="Cort"),   definition = function(copula) {
  # The implementation is in Rcpp.
  bivRho(copula@a,
         copula@b,
         copula@p)
})

#' @describeIn biv_tau-methods Method for the class Cort
setMethod(f = "biv_tau", signature = c(copula="Cort"),   definition = function(copula) {
  # The implmeentation is in Rcpp
  bivTau(copula@a,
         copula@b,
         copula@p)
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

  # The implementation is in Rcpp
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

  # The implementation is in Rcpp
  cpp_result = projectOnTwoDims(a=t(object@a),
                                b=t(object@b),
                                p=object@p,
                                f=object@f,
                                dims=dims,
                                data = object@data[,dims],
                                kern = object@p/object@vols)

  object@dim = 2
  object@data = object@data[,dims]
  object@f = cpp_result$f
  object@p = cpp_result$p
  object@a = cpp_result$a
  object@b = cpp_result$b
  object@vols = cpp_result$vols
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





























