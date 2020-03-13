context("cbCopula tests")
suppressWarnings(library(copula))
library(empCop)
data(LifeCycleSavings)
data(SMI.12)
data(gasoil)
cop <- cbCopula(LifeCycleSavings,m=10)
u=matrix(rep(0,15),ncol=5)
v=matrix(seq(0,1,length.out=15),ncol=5)
w=matrix(rep(1,15),ncol=5)

test_that("zero-row or null data.frame are coerced to indepcopula", {
  expect_error(cbCopula(as.data.frame(NULL)))
  expect_equal(cbCopula(as.data.frame(matrix(0,nrow=0,ncol=5))),indepCopula(5))
})

test_that("data must be provided", {
  expect_error(cbCopula())
})

test_that("Inheritance and methods are there", {
  expect_s4_class(cop,"Copula")
})

test_that("non-dividors m are not allowed", {
  expect_error(cbCopula(LifeCycleSavings,m=3))
  expect_error(cbCopula(LifeCycleSavings,m=7))
  expect_error(cbCopula(LifeCycleSavings,m=11))
  expect_error(cbCopula(LifeCycleSavings,m=49))
})

test_that("dimention of cbCopula is equal dimention of data", {
  expect_equal(dim(cbCopula(LifeCycleSavings)),ncol(LifeCycleSavings))
  expect_equal(dim(cbCopula(SMI.12)),ncol(SMI.12))
  expect_equal(dim(cbCopula(gasoil)),ncol(gasoil))
})

test_that("pCopula values are OK",{
  expect_equal(round(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop),3),c(0,0.024))
  expect_equal(pCopula(matrix(seq(1,1,length.out = 10),nrow=2),cop),c(1,1))
  expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
})

test_that("vCopula did not change",{
  expect_equal(vCopula(u,v,cop)[1],0)
  expect_error(vCopula(v,u,cop))
  expect_equal(vCopula(u,w,cop),rep(1,3))
})

test_that("dim is ok",{
  expect_equal(dim(cop),ncol(LifeCycleSavings))
})

test_that("rCopula output is ok",{
  expect_equivalent(rCopula(0,cbCopula(LifeCycleSavings[,1:2])),rCopula(0,archmCopula("Clayton",0.7,2)))
  expect_is(rCopula(10,cop),"matrix")
  expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

test_that("dCopula returns an error",{
expect_error(dCopula(rep(0.5,dim(cop)),cop))
})


test_that("rCopula output is ok ofr COnvexCombCopula",{
  expect_equivalent(rCopula(0,cop),matrix(ncol=4,nrow=0))
  expect_is(rCopula(10,cop),"matrix")
  expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

test_that("pCopula values are between 0 and 1 with OK bounds.",{
  expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) >= c(0,0)))
  expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) <= c(1,1)))
  expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
  expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop))
})






















