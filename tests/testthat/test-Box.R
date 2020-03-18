context("Box test")
library(cort)

b = Box(a = rep(1/4,4),b = rep(3/4,4))
b2 = Box(rep(0,4),rep(1/2,4))
b3 = Box(rep(0,3),rep(1/2,3))

test_that("initialisation check of Box are ok", {
  testthat::expect_error(Box(c(2,0),c(3,1)))
  testthat::expect_equal(Box(c(0,1/2),c(1,1/2))@volume,0)
  testthat::expect_error(Box(c(1,1,0),c(0,0,1)))
  testthat::expect_error(Box(numeric(),numeric()))
})

test_that("zerobox works", {
  testthat::expect_equal(ZeroBox(3)@dim,3)
  testthat::expect_error(ZeroBox(0))
})

test_that("dimension works", {
  testthat::expect_equal(dim(Box(c(0,0,0,0),c(1,1,1,1))),4)
  testthat::expect_equal(dim(Box(c(0),c(1))),1)
})

test_that("intersection works properly", {
  testthat::expect_equal(intersect(b,b2)@a,rep(1/4,4))
  testthat::expect_equal(intersect(b,b2)@b,rep(1/2,4))
  testthat::expect_error(intersect(b,b3))
})

test_that("simulations and contains works", {
  testthat::expect_equal(nrow(simu_unif(b,100)),100)
  testthat::expect_equal(ncol(simu_unif(b,100)),b@dim)
  testthat::expect_equal(simu_unif(b,0),matrix(0,nrow=0,ncol=b@dim))
  testthat::expect_true(all(contains(b,simu_unif(b,100))))
  testthat::expect_true(contains(b,rep(1/2,b@dim)))
})

test_that("measure_in is OK", {
  testthat::expect_equal(measure_in(b,matrix(0,nrow=7,ncol=4)),rep(0,7))
  testthat::expect_equal(measure_in(b,matrix(0,nrow=0,ncol=4)),numeric())
  testthat::expect_equal(measure_in(b,matrix(1/2,nrow=1,ncol=4)),1/16)
  testthat::expect_equal(measure_in(b,matrix(1,nrow=10,ncol=4)),rep(1,10))
  testthat::expect_error(measure_in(b,rep(0,3)))
  testthat::expect_error(measure_in(b,rep(0,5)))
  testthat::expect_error(measure_in(b,matrix(0,nrow=4,ncol=7)))
  testthat::expect_error(measure_in(b,matrix(0,nrow=0,ncol=7)))
  testthat::expect_error(measure_in(b,matrix(0,nrow=10,ncol=3)))
})


test_that("project is OK", {
  testthat::expect_equal(project(b,c(1,2,3))@dim,3)
  testthat::expect_error(project(b,numeric()))
  testthat::expect_error(project(b,1:10))
  testthat::expect_error(project(b,c(7,8,9)))
  testthat::expect_equal(project(b,c(1,3))@a,b@a[c(1,3)])
})

split = split(b,rep(1/2,2),c(1,2))

test_that("split is OK", {
  testthat::expect_equal(split[[1]]@dim,4)
  testthat::expect_equal(length(split),4)
  testthat::expect_equal(split[[1]]@volume,split[[2]]@volume)
  testthat::expect_error(split(b,rep(1/2,5),c(1,1,2,3,4)))
  testthat::expect_error(split(b,rep(1/2,3),c(1,1,2)))
  testthat::expect_error(split(b,rep(0,3),c(1,2,3)))
  testthat::expect_error(split(b,rep(1,3),c(1,2,3)))
})





















