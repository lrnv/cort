context("Testing the method vCopula")
library(cort)

# define some copulas :
data(LifeCycleSavings)
cop <- cbCopula(LifeCycleSavings[,2:4],m=10)
v=matrix(seq(0,1,length.out=12),ncol=3)
u=matrix(rep(0,12),ncol=3)

testthat::test_that("return 0 for equals inputs", {
  testthat::expect_equal(vCopula(rep(0.5,3),rep(0.5,3),cop),0)
  testthat::expect_equal(vCopula(u,u,cop),rep(0,4))
  testthat::expect_equal(vCopula(v,v,cop),rep(0,4))
})

# test_that("return 1 for complete inputs", {
#   expect_equal(vCopula(rep(0,4),rep(1,4),arch),1)
# })

testthat::test_that("returns error for non-agreing dimentional inputs", {
  testthat::expect_error(vCopula(1,c(1,2),cop))
  testthat::expect_error(vCopula(c(1,2,3),c(1,2),cop))
  testthat::expect_error(vCopula(u,1,cop))
})

