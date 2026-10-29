# Construct an evaluable family specification

A `family_spec` exposes, as pure functions, the pieces koine needs to
compute a pointwise log-likelihood and a posterior-predictive draw for
any lens that does not emit them natively. This low-level constructor is
normally reached through
[`as_family_spec()`](https://max578.github.io/koine/reference/as_family_spec.md),
which adapts a base-R `family` object.

## Usage

``` r
family_spec(
  name,
  link,
  linkinv,
  loglik_pointwise,
  rng,
  support_check,
  aux_params = character(0),
  dlp_deta = NULL,
  d2lp_deta2 = NULL
)
```

## Arguments

- name:

  Family name, e.g. `"gaussian"`, `"binomial"`, `"poisson"`.

- link:

  Link name, e.g. `"identity"`, `"logit"`, `"log"`.

- linkinv:

  Function mapping the linear predictor `eta` to the response-scale mean
  `mu`.

- loglik_pointwise:

  Function `(y, eta, aux)` returning a length-`N` vector of
  `log p(y_i | eta_i, aux)`. The contract koine relies on.

- rng:

  Function `(eta, aux)` returning one posterior-predictive draw
  `y_tilde` (used only when posterior-predictive output is requested).

- support_check:

  Function `(y)` returning `TRUE` if `y` lies in the family support.

- aux_params:

  Character vector naming the auxiliary parameters (e.g. `"sigma"` for
  Gaussian, none for Binomial).

- dlp_deta:

  Optional function `(y, eta, aux)` returning the first derivative of
  the pointwise log-likelihood with respect to the linear predictor (the
  canonical-link score), used by adaptive Gauss-Hermite quadrature to
  locate the random-effect integrand's mode.

- d2lp_deta2:

  Optional function `(y, eta, aux)` returning the second derivative (the
  negative observed information), used for the quadrature's curvature
  scaling.

## Value

An object of class `koine_family_spec`.

## See also

[`as_family_spec()`](https://max578.github.io/koine/reference/as_family_spec.md)

## Examples

``` r
# normally reached through as_family_spec(); the low-level constructor:
fs <- family_spec(
  name = "gaussian", link = "identity",
  linkinv = function(eta) eta,
  loglik_pointwise = function(y, eta, aux) {
    stats::dnorm(y, eta, aux$sigma, log = TRUE)
  },
  rng = function(eta, aux) stats::rnorm(length(eta), eta, aux$sigma),
  support_check = function(y) all(is.finite(y)),
  aux_params = "sigma")
fs$loglik_pointwise(c(0, 1), c(0, 1), list(sigma = 1))
#> [1] -0.9189385 -0.9189385
```
