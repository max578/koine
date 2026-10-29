# Fit the Bayesian ridge lens

Gaussian linear model with a Gaussian (ridge) prior of estimated global
scale via `bayesreg` — mild, near-flat shrinkage anchoring the
low-sparsity end of the shrinkage spectrum. Same resolution rules as
[`fit_horseshoe()`](https://max578.github.io/koine/reference/fit_horseshoe.md).

## Usage

``` r
fit_ridge(
  pass,
  lens_spec = NULL,
  n_iter = 2000L,
  n_burn = 1000L,
  n_chains = 4L
)
```

## Arguments

- pass:

  A `koine_pass` with a Gaussian (identity) family.

- lens_spec:

  The lens's registry record; defaults to the `ridge` record.

- n_iter, n_burn, n_chains:

  Sampler controls.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
d$y <- 1 + 0.8 * d$x1 + rnorm(60)
pass <- koine_pass(y ~ x1 + x2, d, gaussian())
fit_ridge(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
# }
```
