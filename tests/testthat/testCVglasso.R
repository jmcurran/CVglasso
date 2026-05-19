testthat::set.seed(1)

# generate data from a sparse matrix
S <- matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5) {
  for (j in 1:5) {
    S[i, j] <- S[i, j]^abs(i - j)
  }
}

# generate a small n x p matrix with rows drawn from iid N_p(0, S)
Z <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
out <- eigen(S, symmetric = TRUE)
S.sqrt <- out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X <- Z %*% S.sqrt

# common lightweight settings for deterministic and CRAN-light tests
base_args <- list(nlam = 3, K = 2, maxit = 200, trace = "none")

run_and_expect_clean <- function(args) {
  testthat::expect_no_warning(testthat::expect_no_error(do.call(CVglasso, c(args, base_args))))
}

testthat::test_that("CVglasso runs without warnings or errors", {
  run_and_expect_clean(list(X = X, adjmaxit = 100))
  run_and_expect_clean(list(X = X, lam = 0.1, adjmaxit = 100))
  run_and_expect_clean(list(S = S, lam = 0.1, adjmaxit = 100))
})

testthat::test_that("CVglasso parallel options stay lightweight and deterministic", {
  testthat::skip_on_cran()
  run_and_expect_clean(list(X = X, cores = 2, adjmaxit = 100))
  run_and_expect_clean(list(X = X, adjmaxit = 2, cores = 2))
})

testthat::test_that("CVglasso adjmaxit and path options run cleanly", {
  run_and_expect_clean(list(X = X, adjmaxit = 2))
  run_and_expect_clean(list(X = X, path = TRUE, adjmaxit = 100))
  run_and_expect_clean(list(S = S, path = TRUE, adjmaxit = 100))
})
