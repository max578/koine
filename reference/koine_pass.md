# Construct a koine `pass`

Construct a koine `pass`

## Usage

``` r
koine_pass(
  formula,
  data,
  family,
  required_outputs = c("posterior", "log_likelihood"),
  seed = 1L,
  obs_id = NULL,
  resources = list(allowed_tiers = 1L),
  prior = NULL
)
```

## Arguments

- formula:

  A model formula. Population-level terms, optionally with `lme4`-style
  random *intercepts*: a single `(1 | g)`, or the breeder-standard
  multi-environment-trial (MET) structure of two crossed or nested
  random intercepts `(1 | variety) + (1 | trial)`. Richer random-effects
  structure (random slopes, correlated effects, more than two groupings)
  is parsed and flagged unsupported so the battery abstains non-silently
  rather than fitting the wrong model. Whether two groupings are crossed
  or nested is resolved from the data, never from the formula.

- data:

  A data frame holding every column referenced in `formula`.

- family:

  A base-R `family` object or a known family name; wrapped into an
  evaluable
  [`family_spec()`](https://max578.github.io/koine/reference/family_spec.md).

- required_outputs:

  Character subset of
  `c("posterior", "log_likelihood", "posterior_predictive", "sample_stats")`.

- seed:

  Integer master seed; per-lens seeds are derived from it.

- obs_id:

  Optional character vector of observation ids (defaults to
  `obs_0001`-style labels of length `nrow(data)`).

- resources:

  A list with `allowed_tiers` (subset of `1:3`) and optional budgets;
  defaults to tier 1 only.

- prior:

  Optional structured prior table (unused by the DRAM lens, which uses
  its own weakly-informative default; reserved for later lenses).

## Value

A validated object of class `koine_pass`.

## Examples

``` r
d <- data.frame(y = rnorm(20), x = rnorm(20))
koine_pass(y ~ x, d, gaussian())
#> <koine pass>
#>   formula : y ~ x
#>   family  : gaussian / identity
#>   N       : 20
#>   outputs : posterior, log_likelihood
#>   ir_hash : 8ee2601d68fb...
```
