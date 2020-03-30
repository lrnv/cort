#' @include generics.R Box.R WeightedBox.R empiricalCopula.R
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
  root = WeightedBox(Box(rep(0,d),rep(1,d)),data,min_node_size=min_node_size)
  model = .Cort(
    data = data,
    p_value_for_dim_red = p_value_for_dim_red,
    number_max_dim = min(number_max_dim,d),
    min_node_size = min_node_size,
    #leaves=list(root),
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

setMethod(f="fit", signature = c(object="Cort"), definition = function(object,slsqp_options = NULL){

            # Splitting the domain into small boxes :
            if(object@verbose_lvl>0) {cat("Splitting...\n")}

            continue = TRUE
            leaves = list(WeightedBox(Box(rep(0,object@dim),rep(1,object@dim)),object@data,min_node_size=object@min_node_size))

            while(continue){
              are_splittables = purrr::map_lgl(leaves,is_splittable)

              if(object@verbose_lvl>0){cat("\n    ",sum(are_splittables),"leaves to split...")}
              if(object@verbose_lvl>1){cat("\n")}

              if(any(are_splittables)){

                # Split every box :
                leaves = purrr::map(leaves, ~split(.x,
                                  p_val_threshold = object@p_value_for_dim_red,
                                  number_max_dim=object@number_max_dim,
                                  verbose_lvl=object@verbose_lvl-1,
                                  slsqp_options = slsqp_options
                  )) %>%
                  unlist(recursive = FALSE)

                continue = TRUE
              } else {
                continue = FALSE
              }
            }

            # Then compute the weights :
            if(object@verbose_lvl>0) {cat("Enforcing constraints...\n")}

            #  Saving some data :
            object@vols = purrr::map_dbl(leaves,"volume")
            object@f    = purrr::map_dbl(leaves,~nrow(.x@data))/nrow(object@data)
            object@a    = t(sapply(leaves,function(x){x@a})) #  n,d
            object@b    = t(sapply(leaves,function(x){x@b})) # n,d

            # Building constraints :
            n = length(leaves)
            d = ncol(object@data)

            middle_points =  # d,n
            evaluation_points = purrr::map(1:d,~unique((object@b+object@a)[,.x]/2))
            dims = purrr::map2(evaluation_points,1:d,function(x,y){rep(y,length(x))}) %>% unlist()
            F_vec = unlist(evaluation_points)

            lambdas = pmin(pmax((F_vec - t(object@a[,dims]))/(t(object@b[,dims] - object@a[,dims])),0),1)

            # fitting the model :
            model = build_model(P_mat = diag(1/object@vols),
                                q_vec = -object@f/object@vols,
                                A_mat = rbind(lambdas,rep(1,n),diag(n)),
                                l_vec = c(F_vec,1,rep(0,n)),
                                u_vec = c(F_vec,1,rep(Inf,n)),
                                verbose_lvl = object@verbose_lvl-1)

            model$WarmStart(x=object@f)
            rez = model$Solve()

            # saving weights, in the model and in the leaves :
            object@p = pmax(rez$x,0)/sum(pmax(rez$x,0)) # correction for small negative weights.
            leaves = purrr::map2(leaves, object@p,function(l,w){
              l@weight = w
              return(l)
            })
            if(object@verbose_lvl>0){cat("Done !\n")}
            return(object)
          })

build_model <- function(P_mat,q_vec,A_mat,l_vec,u_vec,verbose_lvl=1){

  if (verbose_lvl>0) {
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
  } else {
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

  # remind that pCopula and dCopula generics already transform inputs
  # into matrices...

  if (ncol(u) != dim(copula)) {
    stop("the input value must be coercible to a matrix with dim(copula) columns.")
  }

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
  if (ncol(u) != dim(copula)) {
    stop("the input value must be coercible to a matrix with dim(copula) columns.")
  }

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
      rho[j,i] <- rho[i,j] <- 12*sum(apply((2-b[,c(i,j)] - a[,c(i,j)])/2,1,prod)*p)-3
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

  a <- a[,                rep(1,n_leaves),]
  b <- b[,                rep(1,n_leaves),]
  c <- c[rep(1,n_leaves),                ,]
  d <- d[rep(1,n_leaves),                ,]

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
  cross_B = apply(vapply(1:K,function(i){cross_B[,,dims[i,]]},cross_B[,,dims[1,]]),c(1,2,4),prod) # n,n,K

  # compute boxes measures :
  side_lengths = b - a #n, d
  measures = t(apply(vapply(1:K,function(i){side_lengths[,dims[i,]]},side_lengths[,c(1,2)]),c(1,3),prod)) # K, n

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


  d = ncol(object@data)
  n = nrow(object@a)
  other_n = nrow(other_tree@a)

  kern = object@p/object@vols
  other_kern = other_tree@p/other_tree@vols

  dim(kern) <- c(n,1)
  dim(other_kern) <- c(1,other_n)

  cross_kern = kern[,rep(1,other_n)]*other_kern[rep(1,n),] # n, other_n

  a       = object@a
  b       = object@b
  other_a = other_tree@a
  other_b = other_tree@b

  dim(a) <-       c(n, 1,       d)
  dim(b) <-       c(n, 1,       d)
  dim(other_a) <- c(1, other_n, d)
  dim(other_b) <- c(1, other_n, d)

  a <-             a[,         rep(1,other_n),]
  b <-             b[,         rep(1,other_n),]
  other_a <- other_a[rep(1,n),               ,]
  other_b <- other_b[rep(1,n),               ,]

  side_lengths = pmax(pmin(b,other_b) - pmax(a,other_a),0)
  volumes_of_intersections = apply(side_lengths,1:2,prod) # TAKES A LOT OF TIME

  return(sum(cross_kern*volumes_of_intersections))

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
  colMeans(rng[,rep(1,m)]<=t[rep(1,M),])

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

  leaves = purrr::map(1:new_n,~Box(a = new_a[.x,],b=new_b[.x,]))

  # set things :
  object@data = object@data[,dims]
  object@f = purrr::map_dbl(leaves,function(l){sum(contains(l,object@data))})/nrow(object@data)
  object@p = purrr::map_dbl(leaves,function(l){
      a_complete[dims] <- l@a
      b_complete[dims] <- l@b
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

    graphics::par(mfrow=c(d,d),
      mai = c(0,0,0,0),
      oma = c(3,3,5,3))

  for (i in 1:d){
    for (j in 1:d){
      graphics::plot(c(0,1),c(0,1),type="n",xlab="",ylab="",xaxt='n',yaxt='n')
      if(i != j){
        xx = dd[(dd$dim_x == i)&(dd$dim_y==j),]
        graphics::rect(xx$ymin,xx$xmin,xx$ymax,xx$xmax,col=grDevices::gray(1-xx$col),border=NA,density=NA)
        graphics::points(x@data[,i],x@data[,j],cex=0.7,col="red")
      }
    }
  }

  graphics::par(mfrow = c(1,1))
}





























