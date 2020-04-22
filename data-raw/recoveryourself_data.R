## code to prepare `recoveryourself_data` dataset goes here
set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
x = matrix(runif(1000),500,2)
recoveryourself_data = t(apply(x, 1,function(u){
  if(u[1]< 1/4){
    u[2] = 3/4 + u[2]/4
  } else{ if(u[1]<1/2){
    u[2] = 1/2 + u[2]/4
  } else { if(u[1]<3/4){
    u[2] = u[2]/4
  } else {
    u[2] = 1/4 + u[2]/4
  }}}
  return(u)
}))
usethis::use_data(recoveryourself_data, overwrite = TRUE)
