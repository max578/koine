# Adapt a base-R family to an evaluable `family_spec`

Wraps a base-R / brms `family` object (which exposes `linkinv` but no
pointwise log-likelihood) with the `loglik_pointwise` and `rng`
functions koine needs. Supported families in this build are `gaussian`,
`binomial` (and `bernoulli`, a binomial with one trial), and `poisson`.
An unsupported family raises `error.family_not_evaluable`, never a
silent approximation.

## Usage

``` r
as_family_spec(family, trials = NULL)
```

## Arguments

- family:

  A base-R `family` object (e.g.
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html)), or a single
  string naming one.

- trials:

  Optional length-`N` integer vector of binomial trial counts; defaults
  to 1 (the Bernoulli case).

## Value

A `koine_family_spec`.

## Examples

``` r
fs <- as_family_spec(gaussian())
fs$loglik_pointwise(y = c(0, 1), eta = c(0, 1), aux = list(sigma = 1))
#> [1] -0.9189385 -0.9189385
```
