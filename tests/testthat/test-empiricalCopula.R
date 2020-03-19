context("Testing the class empiricalCopula...")
library(cort)

# we could check that validity is OK.
data("LifeCycleSavings")
cop <- cbCopula(LifeCycleSavings)

cop_wrong <- cop
cop_wrong@pseudo_data <- data.frame("plop")
class(cop_wrong) <- "empiricalCopula"
cop_wrong2 <- cop
class(cop_wrong2) <- "empiricalCopula"
cop_wrong2@pseudo_data <- cop_wrong2@pseudo_data*50

cop_wrong3 <- cop
class(cop_wrong3) <- "empiricalCopula"
cop_wrong3@pseudo_data <- as.data.frame(NULL)

cop_wrong4 <- cop_wrong3
cop_wrong4@pseudo_data <- as.data.frame(matrix(0,nrow=0,ncol=5))

testthat::test_that("validity checks for the empiricalCopula class", {
  testthat::expect_true(validObject(as(cop,"empiricalCopula")))
  testthat::expect_error(supressWarnings(validObject(cop_wrong)))
  testthat::expect_error(validObject(cop_wrong2))
  testthat::expect_error(validObject(cop_wrong3))
  testthat::expect_error(validObject(cop_wrong4))
})
