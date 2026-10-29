# Construct an `fb_fit`

Construct an `fb_fit`

## Usage

``` r
fb_fit(
  draws,
  provenance,
  exactness = c("exact", "approximate", "resampled"),
  sample_stats = NULL,
  log_likelihood = NULL,
  posterior_predictive = NULL,
  diagnostics = NULL,
  status = c("contributed", "abstained"),
  abstention = NULL
)
```

## Arguments

- draws:

  A
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
  of equal-weight posterior draws (F1). Weighted samplers must resample
  to equal weight before constructing the fit and set
  `exactness = "resampled"`.

- provenance:

  A `koine_provenance` record (mandatory, F4).

- exactness:

  `"exact"`, `"approximate"`, or `"resampled"`.

- sample_stats:

  Optional
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
  of HMC-style sample statistics; `NULL` for non-HMC lenses (F2).

- log_likelihood:

  Optional `array[chain, draw, obs_id]`; usually `NULL` at construction
  and filled by
  [`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)
  (F3).

- posterior_predictive:

  Optional predictive draws array.

- diagnostics:

  Optional named list of convergence diagnostics; if `NULL`,
  `rhat`/`ess` are recovered from `draws` via posterior on demand.

- status:

  `"contributed"` or `"abstained"`.

- abstention:

  An `koine_abstention` iff `status == "abstained"`.

## Value

A validated object of class `koine_fb_fit`.

## Examples

``` r
set.seed(1)
draws <- posterior::as_draws_array(
  array(rnorm(8L), dim = c(2L, 2L, 2L),
        dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
prov <- provenance_record(
  method = "demo", method_version = "0", impl_language = "R",
  seed = 1L, derived_seed = 1L, ir_hash = "h",
  maturity_tier = "reference-validated",
  independence_key = list(family = "random_walk", kernel = "demo"),
  role = "opinion", exactness = "exact")
fb_fit(draws = draws, provenance = prov, exactness = "exact",
       status = "contributed")
```
