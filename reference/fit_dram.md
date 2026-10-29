# Fit the DRAM lens

Runs delayed-rejection adaptive Metropolis on the pass's GLM and returns
an [`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md). The
fit carries equal-weight posterior draws and a complete provenance
record; the pointwise log-likelihood is left `NULL` for
[`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)
to compute via the family spec (F3).

## Usage

``` r
fit_dram(pass, lens_spec = NULL, n_iter = 4000L, n_burn = 2000L, n_chains = 4L)
```

## Arguments

- pass:

  A `koine_pass`.

- lens_spec:

  The lens's registry record (for provenance + the derived seed);
  defaults to the registry's `dram` record.

- n_iter:

  Iterations per chain.

- n_burn:

  Burn-in iterations discarded per chain.

- n_chains:

  Number of chains.

## Value

A `koine_fb_fit`.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + 0.8 * d$x + rnorm(40)
pass <- koine_pass(y ~ x, d, gaussian())
fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
# }
```
