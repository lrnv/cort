## code to prepare `impossible_data` dataset goes here
set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
x = matrix(runif(400),200,2)
x = t(apply(x, 1,function(u){
  if(u[1]< 1/3){
    u[2] = 1/2 + u[2]/2
  } else{ if(u[1]<2/3){
    u[2] = u[2]/2
  } else {
    u[2] = 1/2 + u[2]/2
  }}
  return(u)
}))
impossible_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(201)
usethis::use_data(impossible_data, overwrite = TRUE)
