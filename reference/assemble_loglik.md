# Assemble per-lens pointwise log-likelihood matrices (Stage-0 SS6.1)

For each contributed opinion fit, takes the native pointwise
log-likelihood if present, else reconstructs the S x N matrix by
evaluating the family's `loglik_pointwise` at each draw's linear
predictor. Columns are reindexed to the pass's canonical `obs_id` order
(A3).

## Usage

``` r
assemble_loglik(fb_fits, pass)
```

## Arguments

- fb_fits:

  A list of `koine_fb_fit` objects.

- pass:

  The originating `koine_pass`.

## Value

A named list (by lens method) of `list(ll =, chain_id =)`, where `ll` is
an S x N matrix over the canonical observation axis.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40))
pass <- koine_pass(y ~ x, d, gaussian())
fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
ll <- assemble_loglik(list(fit), pass)
dim(ll[[1]]$ll)
#> [1] 600  40
# }
```
