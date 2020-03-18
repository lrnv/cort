context("WeightedBox test")
library(cort)

b = Box(a = rep(1/4,4),b = rep(3/4,4))
wb = WeightedBox(b,matrix(0.5,nrow=10,ncol=4),1,split_dims = 1:4,min_node_size=1)
wb2 = WeightedBox(b,matrix(0.5,nrow=10,ncol=4),1,split_dims = 1:4,min_node_size=11)
wb3 = WeightedBox(b,matrix(0.5,nrow=10,ncol=4),1,split_dims = 1:2,min_node_size=1)

test_that("initialisation check of WeightedBox are ok", {
  testthat::expect_error(WeightedBox(b,0,1,split_dims = 1:4,min_node_size=1))
  testthat::expect_error(WeightedBox(b,matrix(0.5,nrow=4,ncol=3),1,split_dims = 1:4,min_node_size=1))
  testthat::expect_error(WeightedBox(b,matrix(0.5,nrow=7,ncol=4),1,split_dims = 1:4,min_node_size=-3))
  testthat::expect_error(WeightedBox(b,matrix(0,nrow=2,ncol=4),1,split_dims = 1:4,min_node_size=-3))
  testthat::expect_error(WeightedBox(b,matrix(0.5,nrow=4,ncol=3),1,split_dims = 5,min_node_size=-3))
})


test_that("is_splittable works properly", {
  testthat::expect_true(!is_splittable(WeightedBox(b,matrix(0.5,nrow=3,ncol=4),1,split_dims = numeric(),min_node_size=1)))
  testthat::expect_true(!is_splittable(WeightedBox(b,matrix(0.5,nrow=3,ncol=4),1,split_dims = 1:4,min_node_size=6)))
  testthat::expect_true(is_splittable(WeightedBox(b,matrix(0.5,nrow=3,ncol=4),1,split_dims = 1:4,min_node_size=2)))
  testthat::expect_true(!is_splittable(WeightedBox(b,matrix(0.5,nrow=3,ncol=4),1,split_dims = 1,min_node_size=1)))
})

test_that("split behave correctly", {
  testthat::expect_error(split(wb,p_val_threshold=-1,number_max_dim=NULL))
  testthat::expect_error(split(wb,p_val_threshold=-1,number_max_dim=0))
  testthat::expect_true(length(split(wb,p_val_threshold=0,number_max_dim=7))==16)
})
