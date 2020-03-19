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
          definition = function(object,p_val_threshold = 0.75, number_max_dim=object@dim,verbose=FALSE){
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
                                             verbose=verbose)


              bp = optimizer$bp # for the moment we split in the middle.
              p_values = optimizer$p_values

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


                # Compute new boxes and assign the data to each box :
                new_boxes = split(as(object,"Box"),bp,object@split_dims)
                are_the_breakpoint = apply((t(object@data[,object@split_dims]) == bp),2,prod)
                data_without_bp = object@data[!are_the_breakpoint,]

                new_boxes = purrr::map(new_boxes,function(b){
                  dat = object@data[contains(b,data_without_bp,type="loose"),,drop=FALSE]
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

Optmize_breakpoint <- function(data,a=0,b=1,verbose=FALSE){

  # normalise the data # it is transposed
  z = (t(data) - a)/(b-a)
  n = ncol(z)
  d = nrow(z)

  binary_repr = t(sapply(1:(2^d),number2binary,d))

  func = function(bp,t_bin_repr,d,z){
    min = bp*t_bin_repr
    max = bp^(1-t_bin_repr)
    f = purrr::map_dbl(1:2^d, ~mean(apply((min[,.x] <= z)*(z < max[,.x]),2,prod)))
    loss = - sum( (f^2) / apply(max-min,2,prod))
    return(loss)
  }

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
    # arguments to be passed ot the function :
    t_bin_repr = t(binary_repr),
    d = d,
    z = z
  )


  bp = a + optimizer$par*(b-a)
  if(verbose) cat("           breakpoints :",bp,"\n")
  p_values = compute_bootstrapped_p_values(z,bp,binary_repr,n,d)
  if (verbose) cat("           p_values    :",p_values,"\n\n")
  return(list(bp=bp,p_values=p_values))
}

compute_bootstrapped_p_values <- function(z,bp,binary_repr,n,d,N = 199){

  # prerequisites :
  d = nrow(z)
  n = ncol(z)

  min = bp*t(binary_repr)
  max = bp^t(1-binary_repr)

  lambda_l = apply(max-min,2,prod) # 2^d
  lambda_k = vapply(1:d,function(d_rem){lambda_l/((max-min)[d_rem,])},lambda_l) # 2^d * d

  # first, we comput the empirical value of the statistic :
  core = vapply(1:2^d,function(.x){((min[,.x]<=z)*(max[,.x]>z))==1},FUN.VALUE = z) # dims : d, n, 2^d

  f_l = colMeans(apply(core,2:3,prod)) # 2^d
  f_k = vapply(1:d,function(d_rem){ colMeans(apply(core[-d_rem,,,drop=FALSE],2:3,prod))},f_l) # 2^d, d

  statistic <- sum(f_l^2/lambda_l) -2 * colSums(f_k * f_l / lambda_k) # d

  # then we bootstrap it :
  z_rep = vapply(1:N,function(i){z},z) # d, n, N

  z_repeats = vapply(1:d,function(i){
    z = z_rep
    z[i,,] = runif(n*N)
    z
  },z_rep) # d, n, N, D=d

  cores = vapply(1:d,function(d_rem){
    vapply(1:2^d,function(.x){
      ((min[,.x]<=z_repeats[,,,d_rem])*(max[,.x]>z_repeats[,,,d_rem]))==1
    }, FUN.VALUE = z_repeats[,,,d_rem])
  }, FUN.VALUE = array(0,c(d,n,N,2^d))) # d, n, N, 2^d, d

  f_l = colMeans(apply(cores,2:5,prod))
  f_k = vapply(1:d,function(d_rem){
    plop = cores[-d_rem,,,,d_rem,drop=FALSE]
    dim(plop) <- dim(plop)[1:4]
    return(colMeans(apply(plop, 2:4, prod)))
    },array(0.,c(N,2^d))) #(N,2^d,d)

  samples = apply(aperm(f_l^2,c(2,3,1))/lambda_l - 2 * aperm(f_k*f_l,c(2,3,1))/vapply(1:N,function(i){lambda_k},lambda_k),c(2,3),sum)
  p_val = rowMeans(statistic < samples)
  return(p_val)
}






















