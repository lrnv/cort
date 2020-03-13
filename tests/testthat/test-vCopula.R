context("vCopula tests")
library(empCop)

# define some copulas :
arch <- copula::archmCopula('Clayton',0.7,3)
v=matrix(seq(0,1,length.out=12),ncol=3)
u=matrix(rep(0,12),ncol=3)

test_that("return 0 for equals inputs", {
  expect_equal(vCopula(rep(0.5,3),rep(0.5,3),arch),0)
  expect_equal(vCopula(u,u,arch),rep(0,4))
  expect_equal(vCopula(v,v,arch),rep(0,4))
})

# test_that("return 1 for complete inputs", {
#   expect_equal(vCopula(rep(0,4),rep(1,4),arch),1)
# })

test_that("returns error for non-agreing dimentional inputs", {
  expect_error(vCopula(1,c(1,2),arch))
  expect_error(vCopula(c(1,2,3),c(1,2),arch))
  expect_error(vCopula(u,1,arch))
})

