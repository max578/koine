# Exact marginal log-likelihood of a crossed/nested Gaussian MET model (Stage-3)

Computes `log N(y | eta_fixed, V)` with the structured dense marginal
covariance `V = sigma^2 I + sum_r tau_r^2 Z_r Z_r'`, via a Cholesky of
`V`. This is the exact marginal of a Gaussian linear mixed model with
one or more crossed (or nested) random intercepts — the random effects
integrated out, not sampled.

## Usage

``` r
.re_marginal_gaussian_met(y, eta_fixed, factors, taus, sigma)
```

## Arguments

- y:

  Response vector (length `N`).

- eta_fixed:

  Fixed-part linear predictor (length `N`).

- factors:

  The per-grouping list from `.build_re_design()` (each carries the
  integer group index and the group count).

- taus:

  Vector of random-intercept SDs, one per grouping (`>= 0`).

- sigma:

  Residual SD (`> 0`).

## Value

A scalar log marginal likelihood (`-Inf` if `V` is not positive-definite
under the proposed variance components).
