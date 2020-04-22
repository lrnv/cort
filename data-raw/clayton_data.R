## code to prepare `clayton_data` dataset goes here
psi <- function(t,alpha) (1 + sign(alpha)*t) ^ (-1/alpha) # generator
rClayton <- function(n,dim,alpha){
  val <- matrix(runif(n * dim), nrow = n)
  gam <- rgamma(n, shape = 1/alpha, rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  psi(- log(val) / gam,alpha)
}
set.seed(12,kind = "Mersenne-Twister",normal.kind = "Inversion")
clayton_data <- matrix(nrow=200,ncol=4)
clayton_data[,c(1,4,3)] = rClayton(n=200,dim=3,alpha=7)
clayton_data[,2] = runif(200)
clayton_data[,3] <- 1 - clayton_data[,3]
usethis::use_data(clayton_data, overwrite = TRUE)
