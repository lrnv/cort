#' @include generics.R Box.R
NULL

.WeightedBox = setClass(Class = "WeightedBox", contains=c("Box"),
                slots = c(weight = "numeric",
                          data = "matrix",
                          split_dims = "numeric",
                          min_node_size = "numeric"),
                validity = function(object) {
                  errors <- c()
                  if(ncol(object@data) != object@dim){
                    errors <- c(errors,"data should have the same number of columns than dim")
                  }
                  if(!all(object@split_dims %in% 1:object@dim)){
                    errors <- c(errors, "split_dims should be inside the dimensions.")
                  }
                  if (length(errors) == 0)
                    TRUE else errors
                })

WeightedBox = function(box,data,weight=1,split_dims=1:box@dim,min_node_size = 1) {

  # checks :
  if(ncol(data) != box@dim){ stop ("The data shuold have the same umber of colums as the dimension of the box.")}
  if(!all(contains(box,data))){ stop("The box should contain the data")}
  if(!all(split_dims<=box@dim)) {stop("splitting dimensions should be smaller than the box dimensions")}
  if(weight > 1) { stop("The weight should be smaller than 1")}
  if(length(min_node_size) != 1){ stop("The min_node_size should be an integer")}
  if(length(weight) != 1){ stop("The eight should be a single value")}
  if(min_node_size <=0) {stop("The min_node_size should be positive")}


  if(length(weight) != 1){ stop("You should provide only one weight")}


  class(box) <- "WeightedBox"
  box@weight = weight
  box@split_dims = split_dims
  box@data = data
  box@min_node_size = min_node_size
  validObject(box)
  return(box)
}

setMethod(f = "show",
          signature = c(object = "WeightedBox"),
          definition = function(object){
            cat("WeightedBox:(w = ",object@weight,"split_dims=",object@split_dims,")[",object@a,"]-> [",object@b,"]")
          })

setMethod(f = "is_splittable",
          signature = c(object = "WeightedBox"),
          definition = function(object){
            return((length(object@split_dims)>1) && (nrow(object@data)>object@min_node_size))
          })

setMethod(f="split",
          signature = c(object="WeightedBox"),
          definition = function(object,p_val_threshold = 0.75, number_max_dim=object@dim,verbose_lvl=1){
            if(!is_splittable(object)){return(list(object))}

            # Randomize splitting dimensions :
            if(number_max_dim < 2){stop("Splits cannot be done in dimensions smaller than 2...")}
            if(length(object@split_dims) > number_max_dim){
              random_dims = sample(x =object@split_dims,
                                   size=number_max_dim,
                                   replace=FALSE)
              non_taken_dims = object@split_dims[!(object@split_dims %in% random_dims)]
              object@split_dims = random_dims
            } else{
              non_taken_dims = numeric()
            }



            if (is_splittable(object)){ #if it is splittable"
              optimizer = Optmize_breakpoint(object@data[,object@split_dims],
                                             object@a[object@split_dims],
                                             object@b[object@split_dims],
                                             verbose_lvl=verbose_lvl-1)


              bp = optimizer$bp # for the moment we split in the middle.
              p_values = optimizer$p_values

              # if any of the breakpoints are touching the boundary, we remove the dimensions :
              touches_boundary = (((bp <= object@a[object@split_dims])+(bp >= object@b[object@split_dims])))>0
              if(any(touches_boundary)){
                if(verbose_lvl>0){print("Be carrefull, we are removing a splitting dimensions because breakpoint touches the boundary. ")}
                p_values[which(touches_boundary)] = p_val_threshold+1
              }




              # try to remove some pslitting dimensions :
              to_be_removed = p_values > p_val_threshold
              if(all(to_be_removed)){
                object@split_dims = non_taken_dims
                return(list(object))
              } else {

                if(any(to_be_removed)){
                  object@split_dims = object@split_dims[!to_be_removed]
                  bp = bp[!to_be_removed]
                }

                if(length(object@split_dims) == 1){
                  object@split_dims = c(object@split_dims,non_taken_dims)
                  return(list(object))
                }


                # Compute new boxes and assign the data to each box :
                new_boxes = split(as(object,"Box"),bp,object@split_dims)
                are_the_breakpoint = apply((t(object@data[,object@split_dims]) == bp),2,prod)
                data_without_bp = object@data[!are_the_breakpoint,]

                new_boxes = purrr::map(new_boxes,function(b){
                  dat = object@data[contains(b,data_without_bp,type="rllc"),,drop=FALSE]
                  WeightedBox(box = b, data = dat,
                              weight =object@weight * nrow(dat)/nrow(object@data),
                              split_dims = c(object@split_dims,non_taken_dims),
                              min_node_size = object@min_node_size)
                })
                return(new_boxes)
              }

            }
            object@split_dims = c(object@split_dims,non_taken_dims)
            return(list(object))
          })

Optmize_breakpoint <- function(data,a=0,b=1,verbose_lvl=0){

  # Compute prerequisites :
  n = nrow(data)
  d = ncol(data)
  z = (t(data) - a)/(b-a)

  binary_repr = sapply(1:(2^d),number2binary,d)
  dim(binary_repr) = c(d,2^d,1)
  binary_repr = aperm(binary_repr,c(1,3,2))
  binary_repr = binary_repr[,rep(1,n),,drop=FALSE]


  z_rep = z
  dim(z_rep) = c(d,n,1)
  z_rep = z_rep[,,rep(1,2^d),drop=FALSE]

  # The loss function :
  func = function(bp,binary_repr,d,z_rep){
    min = bp*binary_repr
    max = bp^(1-binary_repr)
    f = colMeans(colSums((min <= z_rep)&(z_rep<max))==d)
    xxx = (f^2) / apply(max[,1,]-min[,1,],2,prod)
    loss = - sum( xxx[!is.na(xxx)] ) # permet de ne pas avoir de bug quand il y a des boites de taille 0.
    return(loss)
  }

  # Optimize it :
  optimizer = nloptr::cobyla( # slsqp cobyla, ...
    x0 = rowMeans(z),
    fn = func,
    lower = rep(0,d),
    upper = rep(1,d),
    nl.info=FALSE,
    # control=list(stopval = -Inf, # stop minimization at this value
    #              xtol_rel = 1e-6, # stop on small optimization step
    #              maxeval = 1000, # stop on this many function evaluations
    #              ftol_rel = 0.0, # stop on change times function value
    #              ftol_abs = 0.0, # stop on small change of function value
    #              check_derivatives = FALSE),
    # arguments to be passed to the function :
    binary_repr = binary_repr,
    d = d,
    z_rep = z_rep
  )

  # Return the breakpoint (re-normalized to the box) :
  bp = a + optimizer$par*(b-a)
  if(verbose_lvl>0) cat("           breakpoints :",bp,"\n")

  # Compute the p_values :
  p_values = compute_bootstrapped_p_values(z,bp)
  if (verbose_lvl>0) cat("           p_values    :",p_values,"\n\n")
  return(list(bp=bp,p_values=p_values))
}

compute_bootstrapped_p_values <- function(z,bp,N = 499){

  # prerequisites :
  d = nrow(z)
  n = ncol(z)

  binary_repr = t(sapply(1:(2^d),number2binary,d))
  min = bp*t(binary_repr)
  max = bp^t(1-binary_repr)

  lambda_l = apply(max-min,2,prod) # 2^d
  lambda_k = vapply(1:d,function(d_rem){lambda_l/((max-min)[d_rem,])},lambda_l) # 2^d * d

  # Compute empirical value of the statistic :
  core = vapply(1:2^d,function(.x){((min[,.x]<=z)*(max[,.x]>z))==1},FUN.VALUE = z) # dims : d, n, 2^d

  f_l = colMeans(colSums(core,dims=1)==d) # 2^d
  f_k = vapply(1:d,function(d_rem){ colMeans(colSums(core[-d_rem,,,drop=FALSE],dims=1)==d-1)},f_l) # 2^d, d

  statistic <- sum(f_l^2/lambda_l) -2 * colSums(f_k * f_l / lambda_k) # d

  # NOW BOOTRAP THE SAME THING by first augmenting the dimensionality of z, min, max. The issue is that it eats a bunch of memory an calls GC a lot.
  # maybe this should go in C.
  index = (diag(d)==1)
  dim(index) = c(d,d,1,1,1)
  index = aperm(index,c(1,3,4,5,2))

  dim(z) = c(d,n,1,1,1)
  z = z[,,rep(1,N),,rep(1,d),drop=FALSE]
  z[index[,rep(1,n),rep(1,N),,]] <- runif(n*N*d) # the bootstrap is here.
  z = z[,,,rep(1,2^d),] # TAKES TIME 5

  dim(min) = c(d,2^d,1,1,1)
  min = aperm(min,c(1,3,4,2,5))
  min = min[,rep(1,n),rep(1,N),,rep(1,d)] # TAKES TIME 5

  dim(max) = c(d,2^d,1,1,1)
  max = aperm(max,c(1,3,4,2,5))
  max = max[,rep(1,n),rep(1,N),,rep(1,d)] # TAKES TIME 5

  core_boot = (min <= z) & (max > z) # d, n, N, 2^d, d # TAKES TIME 5

  # Compute bootstrapped_statistic :
  f_l_boot = colMeans(colSums(core_boot,dims=1) == d) # TAKES TIME 1

  f_k_boot = vapply(1:d,function(d_rem){
    colMeans(colSums(core_boot[-d_rem,,,,d_rem,drop=FALSE],dims=1)==(d-1)) # TAKES TIME 5
  },array(0.,c(N,2^d))) #(N,2^d,d)

  # Summarise :
  bootstrap_samples = apply(aperm(f_l_boot^2,c(2,3,1))/lambda_l - 2 * aperm(f_k_boot*f_l_boot,c(2,3,1))/vapply(1:N,function(i){lambda_k},lambda_k),c(2,3),sum) # TAKES TIME 1
  p_val = rowMeans(statistic < bootstrap_samples)
  return(p_val)
}





compute_bootstrapped_p_values_with_rray <- function(z,bp,N = 499){

  # prerequisites :
  d = nrow(z)
  n = ncol(z)

  binary_repr = t(sapply(1:(2^d),number2binary,d))
  min = bp*t(binary_repr)
  max = bp^t(1-binary_repr)

  lambda_l = apply(max-min,2,prod) # 2^d
  lambda_k = vapply(1:d,function(d_rem){lambda_l/((max-min)[d_rem,])},lambda_l) # 2^d * d

  # Compute empirical value of the statistic :
  core = vapply(1:2^d,function(.x){((min[,.x]<=z)*(max[,.x]>z))==1},FUN.VALUE = z) # dims : d, n, 2^d

  f_l = colMeans(colSums(core,dims=1)==d) # 2^d
  f_k = vapply(1:d,function(d_rem){ colMeans(colSums(core[-d_rem,,,drop=FALSE],dims=1)==d-1)},f_l) # 2^d, d

  statistic <- sum(f_l^2/lambda_l) -2 * colSums(f_k * f_l / lambda_k) # d

  # NOW BOOTRAP THE SAME THING by first augmenting the dimensionality of z, min, max. The issue is that it eats a bunch of memory an calls GC a lot.
  # maybe this should go in C.
  index = (diag(d)==1)
  dim(index) = c(d,d,1,1,1)
  index = aperm(index,c(1,3,4,5,2))

  dim(z) = c(d,n,1,1,1)
  z = z[,,rep(1,N),,rep(1,d),drop=FALSE]
  z[index[,rep(1,n),rep(1,N),,]] <- runif(n*N*d) # the bootstrap is here.
  z = z[,,,rep(1,2^d),] # TAKES TIME 5

  dim(min) = c(d,2^d,1,1,1)
  min = aperm(min,c(1,3,4,2,5))
  min = min[,rep(1,n),rep(1,N),,rep(1,d)] # TAKES TIME 5

  dim(max) = c(d,2^d,1,1,1)
  max = aperm(max,c(1,3,4,2,5))
  max = max[,rep(1,n),rep(1,N),,rep(1,d)] # TAKES TIME 5

  core_boot = (min <= z) & (max > z) # d, n, N, 2^d, d # TAKES TIME 5

  # Compute bootstrapped_statistic :
  f_l_boot = colMeans(colSums(core_boot,dims=1) == d) # TAKES TIME 1

  f_k_boot = vapply(1:d,function(d_rem){
    colMeans(colSums(core_boot[-d_rem,,,,d_rem,drop=FALSE],dims=1)==(d-1)) # TAKES TIME 5
  },array(0.,c(N,2^d))) #(N,2^d,d)

  # Summarise :
  bootstrap_samples = apply(aperm(f_l_boot^2,c(2,3,1))/lambda_l - 2 * aperm(f_k_boot*f_l_boot,c(2,3,1))/vapply(1:N,function(i){lambda_k},lambda_k),c(2,3),sum) # TAKES TIME 1
  p_val = rowMeans(statistic < bootstrap_samples)
  return(p_val)
}
















