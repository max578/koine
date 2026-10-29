# Per-group random-effect marginal log-likelihood (Stage-2)

Computes
`sum_g log int_b [prod_{i in g} p(y_i|mu_i(b))] N(b|0,tau^2) db`
(`mu_i(b) = linkinv(eta_fixed_i + b)`) — the marginal likelihood the RE
lenses target when they sample `(beta, sigma, tau)` without ever
sampling the random effects. Exact for the identity-link Gaussian LMM (a
multivariate normal per group, via the determinant / Sherman-Morrison
identities); Gauss-Hermite quadrature otherwise.

## Usage

``` r
.re_group_marginal_loglik(
  family_spec,
  y,
  eta_fixed,
  group,
  tau,
  aux,
  gh = NULL
)
```

## Arguments

- family_spec:

  A `koine_family_spec`.

- y:

  Response vector (length `N`).

- eta_fixed:

  Fixed-part linear predictor (length `N`).

- group:

  Integer group index per observation (length `N`).

- tau:

  Random-intercept standard deviation (scalar, `>= 0`).

- aux:

  Named list of auxiliary parameters (e.g. `list(sigma = .)`).

- gh:

  A Gauss-Hermite rule from `.gh_rule()`; computed if `NULL`.

## Value

A scalar log marginal likelihood.
