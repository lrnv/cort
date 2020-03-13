context("ConvexCombCopula tests")
suppressWarnings(library(copula))
library(empCop)

## ----fig.cap="constructing a copula"-----------------------
copulas <- list(
  copula::archmCopula('gumbel',3,dim=2),
  copula::archmCopula('clayton',-1,dim=2)
)
alpha <- c(1,4)

cop <- ConvexCombCopula(copulas,alpha)

test_that("inputed argument for copulas is handeld in the right way.", {
  expect_error(ConvexCombCopula(copulas = "abc"))
  expect_equal(suppressWarnings(ConvexCombCopula(copulas[[1]])),copulas[[1]])
  expect_equal(suppressWarnings(ConvexCombCopula(copulas[1])),copulas[[1]])
  expect_warning(ConvexCombCopula(copulas[[1]]))
  expect_warning(ConvexCombCopula(copulas[1]))
})

test_that("ConvexCombCopula class herited properly", {
  expect_is(cop,"ConvexCombCopula")
  expect_is(cop,"Copula")
})

test_that("method dim for ConvexCombCopula class works properly", {
  expect_equal(dim(cop),dim(copulas[[1]]))
})

test_that("rCopula output is ok for COnvexCombCopula",{
  expect_equivalent(rCopula(0,cop),matrix(ncol=dim(cop),nrow=0))
  expect_is(rCopula(10,cop),"matrix")
  expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

test_that("pCopula values are between 0 and 1 with OK bounds.",{
  expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 4),nrow=2),cop) >= c(0,0)))
  expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 4),nrow=2),cop) <= c(1,1)))
  expect_equal(pCopula(matrix(seq(0,0,length.out = 4),nrow=2),cop),c(0,0))
  expect_error(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop))
})