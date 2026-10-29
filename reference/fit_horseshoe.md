# Fit the horseshoe (sparsity) lens

Gaussian linear model with a horseshoe prior (Carvalho, Polson & Scott
2010) via `bayesreg`. Resolves only for a continuous (Gaussian,
identity) response with at least one predictor; otherwise a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md).

## Usage

``` r
fit_horseshoe(
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

  The lens's registry record; defaults to the `horseshoe` record.

- n_iter, n_burn, n_chains:

  Sampler controls.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## References

Carvalho, Polson & Scott (2010). The horseshoe estimator for sparse
signals. *Biometrika* 97(2):465-480.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
d$y <- 1 + 0.8 * d$x1 + rnorm(60)
pass <- koine_pass(y ~ x1 + x2, d, gaussian())
fit_horseshoe(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
# }
```
