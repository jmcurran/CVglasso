## ----r setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = FALSE, tidy = TRUE)

## ----r, message = FALSE----------------------------------------------------
library(CVglasso)
library(microbenchmark)

#  generate data from tri-diagonal (sparse) matrix
# compute covariance matrix (can confirm inverse is tri-diagonal)
S = matrix(0, nrow = 100, ncol = 100)

for (i in 1:100){
  for (j in 1:100){
    S[i, j] = 0.7^(abs(i - j))
  }
}

# generate 1000 x 100 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(1000*100), nrow = 1000, ncol = 100)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance matrix
sample = (nrow(X) - 1)/nrow(X)*cov(X)


## ----r, message = FALSE----------------------------------------------------

# benchmark CVglasso - defaults
microbenchmark(CVglasso(S = sample, lam = 0.1, trace = "none"))


## ----r, message = FALSE----------------------------------------------------

# benchmark CVglasso - tolerance 1e-6
microbenchmark(CVglasso(S = sample, lam = 0.1, tol = 1e-6, trace = "none"))


## ----r, message = FALSE----------------------------------------------------

# benchmark CVglasso CV - default parameter grid
microbenchmark(CVglasso(X, trace = "none"), times = 5)


## ----r, message = FALSE----------------------------------------------------

# benchmark CVglasso parallel CV
microbenchmark(CVglasso(X, cores = 2, trace = "none"), times = 5)
