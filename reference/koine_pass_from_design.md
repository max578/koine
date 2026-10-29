# Construct a koine `pass` from design arrays (the IR-targetable builder)

The programmatic, re-entrant, deparse-free entry point a parent
triangulation panel drives (flexyBayes backend-contract clause 1): build
a pass directly from a named design matrix, a response, a family, and an
optional grouping vector — no model formula, no
[`match.call()`](https://rdrr.io/r/base/match.call.html) capture, no
session-global state. Two passes can be built and synthesised in the
same process without interference. The stored design is used verbatim
downstream, so non-syntactic canonical column names (e.g. `(Intercept)`)
pass through untouched.

## Usage

``` r
koine_pass_from_design(
  x,
  y,
  family,
  group = NULL,
  groups = NULL,
  coef_names = colnames(x),
  response_name = "y",
  grouping_name = "g",
  required_outputs = c("posterior", "log_likelihood"),
  seed = 1L,
  obs_id = NULL,
  resources = list(allowed_tiers = 1L),
  prior = NULL
)
```

## Arguments

- x:

  A numeric design matrix with named columns (the fixed-effects design;
  include an `(Intercept)` column if wanted). Column names become the
  coefficient names (`b_<col>` in the draws).

- y:

  The response vector (length `nrow(x)`).

- family:

  A base-R `family` object or a known family name.

- group:

  Optional grouping vector (length `nrow(x)`) for a single random
  intercept `(1 | group)`; `NULL` for a population model. For the MET
  structure of two crossed/nested random intercepts pass `groups`
  instead.

- groups:

  Optional named list of two grouping vectors (each length `nrow(x)`)
  for two crossed/nested random intercepts — the breeder-standard MET
  structure, e.g. `list(variety = ., trial = .)`. The list names become
  the grouping-factor names; mutually exclusive with `group`.

- coef_names:

  Optional column names for `x` (defaults to `colnames(x)`).

- response_name:

  Name recorded for the response (default `"y"`).

- grouping_name:

  Name recorded for the (single) grouping factor when `group` is
  supplied (default `"g"`); ignored when `groups` is used.

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

## See also

[`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md),
[`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md)

## Examples

``` r
set.seed(1)
x <- cbind(`(Intercept)` = 1, z = rnorm(40))
y <- 1 + 0.5 * x[, "z"] + rnorm(40)
pass <- koine_pass_from_design(x, y, gaussian())
pass$has_re             # FALSE -- a population model, built without a formula
#> [1] FALSE
```
