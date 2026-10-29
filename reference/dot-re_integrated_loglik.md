# Per-observation integrated random-effect predictive (Stage-2)

Computes
`log integral p(y_i | linkinv(eta_fixed_i + b)) N(b | 0, tau^2) db` for
every observation — the predictive density for a *new* observation with
its random effect marginalised out. Analytic for the identity-link
Gaussian LMM; Gauss-Hermite quadrature otherwise. This is the quantity
every RE lens is stacked on, so the LOO weights are comparable across
lenses (the same held-out target).

## Usage

``` r
.re_integrated_loglik(family_spec, y, eta_fixed, tau, aux, gh = NULL)
```

## Arguments

- family_spec:

  A `koine_family_spec`.

- y:

  Response vector (length `N`).

- eta_fixed:

  Fixed-part linear predictor `x_i' beta` (length `N`).

- tau:

  Random-intercept standard deviation (scalar, `>= 0`).

- aux:

  Named list of auxiliary parameters (e.g. `list(sigma = .)`).

- gh:

  A Gauss-Hermite rule from `.gh_rule()`; computed if `NULL`.

## Value

A length-`N` vector of integrated pointwise log-likelihoods.
