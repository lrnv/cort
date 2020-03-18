context("Tree test")
library(cort)

model = Cort(LifeCycleSavings)
fitted_model = fit(model)

test_that("initialisation check of Tree are ok", {
  testthat::expect_error(model = Cort(LifeCycleSavings,pseudo_data=TRUE))
})
