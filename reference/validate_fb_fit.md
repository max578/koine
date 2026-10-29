# Validate an `fb_fit` against the field rules

Enforces: a contributed fit carries draws + a complete provenance tuple
(F4); `exactness` agrees between the fit and its provenance; an
abstained fit carries an `koine_abstention` and no draws;
`sample_stats`/`log_likelihood` are either `NULL` or arrays (no
fabricated values, I1).

## Usage

``` r
validate_fb_fit(x)
```

## Arguments

- x:

  A `koine_fb_fit`.

## Value

`x`, invisibly, if valid; otherwise stops.

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
validate_fb_fit(fb_fit(draws, prov, "exact", status = "contributed"))
```
