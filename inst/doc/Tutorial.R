## ----r setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = FALSE, tidy = TRUE)

## ----r, message = FALSE, echo = TRUE---------------------------------------
library(CVglasso)

#  generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# snap shot of data
head(X)


## ----r, message = FALSE, echo = TRUE---------------------------------------

# print oracle covariance matrix
S

# print inverse covariance matrix (omega)
round(qr.solve(S), 5)


## ----r, message = FALSE, echo = TRUE---------------------------------------

# print inverse of sample precision matrix (perhaps a bad estimate)
round(qr.solve(cov(X)*(nrow(X) - 1)/nrow(X)), 5)


## ----r, message = FALSE, echo = TRUE---------------------------------------

# cross validation for lam
CVglasso(X, trace = "none")


## ----r, message = FALSE, echo = TRUE---------------------------------------

# produce CV heat map
CV = CVglasso(X, trace = "none")
plot(CV, type = "heatmap")


## ----r, message = FALSE, echo = TRUE---------------------------------------

# produce line graph for CV errors
plot(CV, type = "line")


## ----r, message = FALSE, echo = TRUE---------------------------------------

# AIC
plot(CVglasso(X, crit.cv = "AIC", trace = "none"))

# BIC
plot(CVglasso(X, crit.cv = "BIC", trace = "none"))


## ----r, message = FALSE, echo = TRUE---------------------------------------

# keep all estimates using path
CV = CVglasso(X, path = TRUE, trace = "none")

# print only first three objects
CV$Path[,,1:3]


## ----r, message = FALSE, echo = TRUE, eval = FALSE-------------------------

# reduce number of lam to 5
CV = CVglasso(X, nlam = 5)


## ----r, message = FALSE, echo = TRUE, eval = FALSE-------------------------

# reduce number of folds to 3
CV = CVglasso(X, K = 3)


## ----r, message = FALSE, echo = TRUE, eval = FALSE-------------------------

# relax convergence criteria
CV = CVglasso(X, tol = 1e-3)


## ----r, message = FALSE, echo = TRUE, eval = FALSE-------------------------

# adjust maximum number of iterations
CV = CVglasso(X, maxit = 1e3)


## ----r, message = FALSE, echo = TRUE, eval = FALSE-------------------------

# parallel CV
CV = CVglasso(X, cores = 3)
