## code to prepare `funcdep_data` dataset goes here
# getting parameters :
set.seed(seed = 12,kind = "Mersenne-Twister",normal.kind = "Inversion")
x = matrix(runif(1500),500,3)
x[,2] = sin(2*pi*x[,1])-x[,2]/pi
x[,3] = (x[,3]*(x[,1]<1/4)/2 - sin(pi**(x[,1]))*(x[,1]>1/4))*(1+x[,3]/(pi^2))
funcdep_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(501)
usethis::use_data(funcdep_data, overwrite = TRUE)
