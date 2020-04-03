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

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


toString.data.frame = function (object, ..., digits=NULL, quote=FALSE, right=TRUE, row.names=TRUE) {

  # taken and adapted from there : https://stackoverflow.com/a/45541857/8425270

  nRows = length(row.names(object));
  # get text-formatted version of the data.frame
  m = as.matrix(format.data.frame(object, digits=digits, na.encode=FALSE));
  # define row-names (if required)
  if (isTRUE(row.names)) {
    rowNames = dimnames(object)[[1]];
    # add empty header (used with column headers)
    rowNames = c("", rowNames);
  }
  # add column headers
  m = rbind(dimnames(m)[[2]], m);
  # add row headers
  m = cbind(rowNames, m);
  # max-length per-column
  maxLen = apply(apply(m, c(1,2), nchar), 2, max, na.rm=TRUE);

  # add right padding
  ##  t is needed because "If each call to FUN returns a vector of length n, then apply returns an array of dimension c(n, dim(X)[MARGIN])"
  m = t(apply(m, 1, str_pad, width=maxLen, side="right"));
  m = t(apply(m, 1, str_pad, width=maxLen+3, side="left"));
  # merge columns
  m = apply(m, 1, paste, collapse="");
  # merge rows (and return)
  return(paste(m, collapse="\n"));
}

str_pad <- function(string,width,side){ # retro-engeniered from stringr::str_pad to remove the dependence...
  mapply(function(str,width){
    n = nchar(str)
    if(side == "left")  return(paste0(paste0(rep(" ",max(width-n,0)),collapse=""),str))
    if(side == "right") return(paste0(str,paste0(rep(" ",max(width-n,0)),collapse="")))
  },string,width,USE.NAMES = FALSE)
}

