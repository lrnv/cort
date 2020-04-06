# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cortMonteCarlo <- function(z, min, max, N) {
    .Call(`_cort_cortMonteCarlo`, z, min, max, N)
}

quadProd <- function(a, b, kern, other_a, other_b, other_kern) {
    .Call(`_cort_quadProd`, a, b, kern, other_a, other_b, other_kern)
}

normMatrix <- function(as, bs, kernels) {
    .Call(`_cort_normMatrix`, as, bs, kernels)
}

lossFunc <- function(bp, bin_repr, z) {
    .Call(`_cort_lossFunc`, bp, bin_repr, z)
}
