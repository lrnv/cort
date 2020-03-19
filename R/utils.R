#### A function to normalise the data that is imputed to methods of the package
normalise_data <- function(u,dim=NULL){

  if (!is.matrix(u))
    u <- rbind(u, deparse.level = 0L)
  if(!is.null(dim)){
    if(ncol(u) != dim){
      stop("You should provide a matrix with the right number of columns")
    }
  }
  ## here as well, 'outside' and 'on-boundary' are equivalent:
  u[] <- pmax(0, pmin(1, u))
  return (u)
}

# convert integrer to binary format
number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  binary_vector[-(1:(length(binary_vector) - noBits))]
}

# A `sample` function more efficient (cf ?sample)
resample <- function(x, ...) x[sample.int(length(x), ...)]


boxes_from_points <- function(points,m){
  t(floor(t(points)*m)/m)
}






