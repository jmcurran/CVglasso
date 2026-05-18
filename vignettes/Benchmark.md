---
title: "Benchmarks"
author: "Matt Galloway"
#date: "2026-05-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Benchmarks}
  %\VignetteEngine{knitr::knitr}
  %\usepackage[UTF-8]{inputenc}
---



This is a short effort to give users an idea of how long the functions take to process. The benchmarks were performed using the default R install on [Travis CI](https://travis-ci.org/).

We will be estimating a tri-diagonal precision matrix with dimension $p = 100$:

<br>\vspace{0.5cm}

``` r
library(CVglasso)
library(microbenchmark)

# generate data from tri-diagonal (sparse) matrix compute covariance matrix
# (can confirm inverse is tri-diagonal)
S = matrix(0, nrow = 100, ncol = 100)

for (i in 1:100) {
    for (j in 1:100) {
        S[i, j] = 0.7^(abs(i - j))
    }
}

# generate 1000 x 100 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(1000 * 100), nrow = 1000, ncol = 100)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance matrix
sample = (nrow(X) - 1)/nrow(X) * cov(X)
```
<br>\vspace{0.5cm}

 - Default convergence tolerance with specified tuning parameter (no cross validation):

<br>\vspace{0.5cm}

``` r
# benchmark CVglasso - defaults
microbenchmark(CVglasso(S = sample, lam = 0.1, trace = "none"))
```

```
## Warning in microbenchmark(CVglasso(S = sample, lam = 0.1, trace = "none")):
## less accurate nanosecond times to avoid potential integer overflows
```

```
## Unit: milliseconds
##                                             expr      min       lq     mean
##  CVglasso(S = sample, lam = 0.1, trace = "none") 22.29949 22.50338 24.27309
##    median       uq      max neval
##  23.05483 23.32459 72.18993   100
```
<br>\vspace{0.5cm}

 - Stricter convergence tolerance with specified tuning parameter (no cross validation):

<br>\vspace{0.5cm}

``` r
# benchmark CVglasso - tolerance 1e-6
microbenchmark(CVglasso(S = sample, lam = 0.1, tol = 1e-06, trace = "none"))
```

```
## Unit: milliseconds
##                                                          expr      min       lq
##  CVglasso(S = sample, lam = 0.1, tol = 1e-06, trace = "none") 38.73819 38.81058
##      mean   median       uq      max neval
##  39.51906 38.96197 39.90052 45.30217   100
```
<br>\vspace{0.5cm}

 - Default convergence tolerance with cross validation for `lam`:

<br>\vspace{0.5cm}

``` r
# benchmark CVglasso CV - default parameter grid
microbenchmark(CVglasso(X, trace = "none"), times = 5)
```

```
## Unit: milliseconds
##                         expr      min      lq     mean   median       uq
##  CVglasso(X, trace = "none") 810.9902 833.287 854.8546 851.9119 864.9527
##       max neval
##  913.1314     5
```
<br>\vspace{0.5cm}

 - Parallel (`cores = 2`) cross validation:

<br>\vspace{0.5cm}


