library(cort)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("should be fine", {
  expect_equal(4,4)
})

test_that("pipe works correctly", {
  expect_equal(c(1,2) %>% sum,3)
})
