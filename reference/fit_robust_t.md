# Fit the robust-t lens (Student-t linear model)

Fits `y ~ Student-t(nu, X beta, sigma)` by delayed-rejection adaptive
Metropolis (the DRAM kernel) and returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md)
carrying its own Student-t pointwise log-likelihood. Resolves only for a
continuous (Gaussian-declared) response; otherwise returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md).

## Usage

``` r
fit_robust_t(
  pass,
  lens_spec = NULL,
  nu = 4,
  n_iter = 4000L,
  n_burn = 2000L,
  n_chains = 4L
)
```

## Arguments

- pass:

  A `koine_pass` with a Gaussian (identity) family.

- lens_spec:

  The lens's registry record; defaults to the `robust_t` record.

- nu:

  Degrees of freedom of the Student-t likelihood (fixed; lower = of
  heavier tails). Default 4.

- n_iter, n_burn, n_chains:

  Sampler controls.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x = rnorm(50))
d$y <- 1 + 0.8 * d$x + stats::rt(50, df = 4) * 0.6
pass <- koine_pass(y ~ x, d, gaussian())
fit_robust_t(pass, nu = 4, n_iter = 800L, n_burn = 400L, n_chains = 2L)
# }
```
