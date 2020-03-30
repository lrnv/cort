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
          definition = function(object,p_val_threshold = 0.75, number_max_dim=object@dim,verbose_lvl=1,slsqp_options=NULL){

            if(!is_splittable(object)){return(list(object))}

            # prepare the verbosity d.f :
            verb_df = data.frame(min = object@a, max = object@b)
            verb_df$bp      = rep(NaN,ncol(object@data))
            verb_df$p_value = rep(NaN,ncol(object@data))
            verb_df$action  = rep("",ncol(object@data))
            verb_df$reason  = rep("",ncol(object@data))
            row.names(verb_df) = paste0("             ",1:nrow(verb_df))

            n = nrow(object@data)

            if(verbose_lvl>1){cat(paste0("        Leaf with ",n," points.\n"))}

            # Randomize splitting dimensions :
            if(number_max_dim < 2){stop("Splits cannot be done in dimensions smaller than 2...")}
            if(length(object@split_dims) > number_max_dim){
              random_dims       = sample(x =object@split_dims, size=number_max_dim, replace=FALSE)
              non_taken_dims    = object@split_dims[!(object@split_dims %in% random_dims)]
              object@split_dims = random_dims
            } else{
              non_taken_dims = numeric()
            }

            verb_df$action[non_taken_dims] = "Dissmissed"
            verb_df$reason[non_taken_dims] = "Randomly"

            if (is_splittable(object)){

              optimizer = Optmize_breakpoint(object@data[,object@split_dims],
                                             object@a[object@split_dims],
                                             object@b[object@split_dims],
                                             verbose_lvl=verbose_lvl-1,
                                             slsqp_options = slsqp_options)

              verb_df$bp[object@split_dims]      <- bp       <- optimizer$bp
              verb_df$p_value[object@split_dims] <- p_values <- optimizer$p_values

              # if any of the breakpoints are too close to boundary, we remove the dimensions :
              normed_bp      = (bp - object@a[object@split_dims])/(object@b[object@split_dims] - object@a[object@split_dims])
              seuil          = 1/min((n+1)^2,2*n+1,1000)
              close_to_bound = (normed_bp< seuil) + (normed_bp > 1-seuil) > 0
              p_val_too_big  = p_values > p_val_threshold


              verb_df$action[object@split_dims[close_to_bound]] = "Removed"
              verb_df$reason[object@split_dims[close_to_bound]] = "Close to boundary"
              verb_df$action[object@split_dims[p_val_too_big]] = "Removed"
              verb_df$reason[object@split_dims[p_val_too_big]] = "Independence test"

              to_be_removed = p_val_too_big+close_to_bound>0

              if(all(to_be_removed)){
                object@split_dims = non_taken_dims
                result = list(object)
              } else {

                if(any(to_be_removed)){
                  object@split_dims = object@split_dims[!to_be_removed]
                  bp = bp[!to_be_removed]
                }

                if(length(object@split_dims) == 1){

                  verb_df$action[object@split_dims] = "Dissmissed"
                  verb_df$reason[object@split_dims] = "No one-dim split" # NOW the df is full.
                  object@split_dims = c(object@split_dims,non_taken_dims)
                  result = list(object)
                } else {
                  # Compute new boxes and assign the data to each box :
                  verb_df$action[object@split_dims] = "Splitted"

                  new_boxes          = split(as(object,"Box"),bp,object@split_dims)
                  are_the_breakpoint = apply((t(object@data[,object@split_dims]) == bp),2,prod)
                  data_without_bp    = object@data[!are_the_breakpoint,]

                  result = purrr::map(new_boxes,function(b){
                    dat = data_without_bp[contains(b,data_without_bp,type="rllc"),,drop=FALSE]
                    WeightedBox(box = b, data = dat,
                                weight =object@weight * nrow(dat)/n,
                                split_dims = c(object@split_dims,non_taken_dims),
                                min_node_size = object@min_node_size)
                  })
                }
              }
            } else{
              object@split_dims = c(object@split_dims,non_taken_dims)
              result = list(object)
            }

            if(verbose_lvl>2) {
              cat(toString.data.frame(verb_df,digits=8))
              cat("\n")
            }
            if(verbose_lvl>1){cat("\n")}
            return(result) #################################### RETURN
          })

Optmize_breakpoint <- function(data,a=0,b=1,verbose_lvl=0,slsqp_options=NULL){

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

  DEFAULT_SLQP_OPTIONS = list(
    stopval = -Inf,
    xtol_rel = 1e-4,
    maxeval = 100000,
    ftol_rel = 1e-6,
    ftol_abs = 1e-6
  )

  # get setted options :
  if(is.null(slsqp_options)){
    slsqp_options = DEFAULT_SLQP_OPTIONS
  }
  if(!is.null(slsqp_options)){
    if(is.null(slsqp_options$stopval))  slsqp_options$stopval  = DEFAULT_SLQP_OPTIONS$stopval
    if(is.null(slsqp_options$xtol_rel)) slsqp_options$xtol_rel = DEFAULT_SLQP_OPTIONS$xtol_rel
    if(is.null(slsqp_options$maxeval))  slsqp_options$maxeval  = DEFAULT_SLQP_OPTIONS$maxeval
    if(is.null(slsqp_options$ftol_rel)) slsqp_options$ftol_rel = DEFAULT_SLQP_OPTIONS$ftol_rel
    if(is.null(slsqp_options$ftol_abs)) slsqp_options$ftol_abs = DEFAULT_SLQP_OPTIONS$ftol_abs
  }

  optimizer = nloptr::slsqp(
      x0 = rowMeans(z),
      fn = func,
      lower = rep(0,d),
      upper = rep(1,d),
      nl.info=verbose_lvl>2,
      control=slsqp_options,
      binary_repr = binary_repr,
      d = d,
      z_rep = z_rep)

  # Return the breakpoint (re-normalized to the box) :
  bp = a + optimizer$par*(b-a)
  # if(verbose_lvl>0) cat("           breakpoints :",bp,"\n")

  # Compute the p_values :
  p_values = compute_bootstrapped_p_values(z,bp)
  # if (verbose_lvl>0) cat("           p_values    :",p_values,"\n")
  return(list(bp=bp,p_values=p_values))
}

compute_bootstrapped_p_values <- function(z,bp,N = 999){

  # prerequisites :
  d = nrow(z)
  n = ncol(z)

  binary_repr = t(sapply(1:(2^d),number2binary,d))
  min = bp*t(binary_repr)
  max = bp^t(1-binary_repr)

  # calling the cpp function to do that : it's much faster...
  cpp_result = cortMonteCarlo(z,min,max,as.integer(N))
  p_val = rowMeans(cpp_result[1,] <= t(cpp_result[-1,]))

  if(any(is.na(p_val))){
    p_val[is.na(p_val)] = 0
  }

  return(p_val)
}






