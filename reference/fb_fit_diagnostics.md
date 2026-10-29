# Recover convergence diagnostics from an `fb_fit`'s draws

Computes split-Rhat and bulk/tail ESS from the posterior draws via
`posterior::`, engine-agnostically. Never fabricates an HMC-only
statistic (I1, F2): for non-HMC lenses `sample_stats` stays `NULL` and
only the draws-derived diagnostics are returned.

## Usage

``` r
fb_fit_diagnostics(x)
```

## Arguments

- x:

  A contributed `koine_fb_fit`.

## Value

A named list with `rhat`, `ess_bulk`, `ess_tail` (per variable).

## Examples

``` r
set.seed(1)
draws <- posterior::as_draws_array(
  array(rnorm(80L), dim = c(20L, 2L, 2L),
        dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
prov <- provenance_record(
  method = "demo", method_version = "0", impl_language = "R",
  seed = 1L, derived_seed = 1L, ir_hash = "h",
  maturity_tier = "reference-validated",
  independence_key = list(family = "random_walk", kernel = "demo"),
  role = "opinion", exactness = "exact")
fit <- fb_fit(draws = draws, provenance = prov, exactness = "exact",
              status = "contributed")
fb_fit_diagnostics(fit)
#> Warning: The ESS has been capped to avoid unstable estimates.
#> $rhat
#> b_(Intercept)           b_x 
#>     0.9679600     0.9908334 
#> 
#> $ess_bulk
#> b_(Intercept)           b_x 
#>      46.87486      64.08240 
#> 
#> $ess_tail
#> b_(Intercept)           b_x 
#>      49.57265      26.47059 
#> 
```
