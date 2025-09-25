# BiBayMR

Bidirectional Bayesian Mendelian Randomization (BiBayMR) with RcppArmadillo.

## Install

```r
install.packages(c("devtools","Rcpp","RcppArmadillo","roxygen2"))
devtools::document()
Rcpp::compileAttributes()
devtools::install()
```

## Example

```r
set.seed(123)
p <- 50
theta_XY_true <- 0.2; theta_YX_true <- 0.0
gamma <- rnorm(p, 0, 0.05)
delta <- rnorm(p, 0, 0.05)
bX <- gamma + rnorm(p, 0, 0.05)
bY <- theta_XY_true * gamma + delta + rnorm(p, 0, 0.05)
sX <- rep(0.05, p); sY <- rep(0.05, p)

fit <- BiBayMR(bX, sX, bY, sY, n_iter = 2000, burnin = 500, thin = 2, seed = 1)
print(fit)

boot <- BiBayMR_DP(bX, sX, bY, sY, B = 20, n_iter = 1500, burnin = 500, thin = 2, seed = 2)
boot$summary
```
