context("Testing the class Cort...")
library(cort)
model = Cort(LifeCycleSavings[,1:3],verbose=FALSE)

testthat::test_that("initialisation check of Cort are ok", {
  testthat::expect_error(model = Cort(LifeCycleSavings,pseudo_data=TRUE,verbose=FALSE))
})
