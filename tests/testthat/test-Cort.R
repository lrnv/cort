context("Testing the class Cort...")
library(cort)
model = Cort(LifeCycleSavings[,1:3],verbose_lvl = 0)

testthat::test_that("initialisation check of Cort are ok", {
  testthat::expect_error(model = Cort(LifeCycleSavings,pseudo_data=TRUE,verbose_lvl = 0))
})

testthat::test_that("the quandratic norm is coherent with the quad prod", {
  testthat::expect_equal(quad_prod(model,model),quad_norm(model))
})
