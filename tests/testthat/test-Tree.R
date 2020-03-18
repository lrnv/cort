context("Tree test")
library(cort)

test_that("initialisation check of Tree are ok", {
  testthat::expect_error(model = Cort(LifeCycleSavings,pseudo_data=TRUE,verbose=FALSE))
})
