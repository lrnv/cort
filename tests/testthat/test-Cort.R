context("Testing the class Cort...")
library(cort)
model = Cort(LifeCycleSavings,verbose=FALSE)

testthat::test_that("initialisation check of Cort are ok", {
  testthat::expect_error(model = Cort(LifeCycleSavings,pseudo_data=TRUE,verbose=FALSE))
})
