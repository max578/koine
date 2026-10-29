# Can koine contribute a diverse opinion for this model? (pre-flight)

A cheap capability predicate — it resolves the battery but fits nothing
— for a parent panel's backend-availability gate. Returns `TRUE` when at
least two independent opinion lenses would enrol, or a single reason
string from the abstention taxonomy otherwise (`"feature_unsupported"`,
`"no_diverse_battery"`, or `"pass_invalid"` for a malformed model).
Because the PSIS reliability gate and weight-degeneracy check can only
be evaluated after fitting, a `TRUE` here is necessarily optimistic:
koine may still abstain at
[`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)
time. The predicate is a *floor*, not a guarantee.

## Usage

``` r
koine_can_contribute(
  formula,
  data,
  family,
  resources = list(allowed_tiers = 1L),
  registry = koine_stage1_registry()
)
```

## Arguments

- formula, data, family:

  As for
  [`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md).

- resources:

  A resources list (see
  [`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md));
  defaults to tier 1.

- registry:

  A registry list; defaults to
  [`koine_stage1_registry()`](https://max578.github.io/koine/reference/koine_stage1_registry.md).

## Value

`TRUE`, or a single reason-code string.

## See also

[`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md),
[`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)

## Examples

``` r
d <- data.frame(y = rnorm(40), x = rnorm(40))
koine_can_contribute(y ~ x, d, gaussian())
#> [1] TRUE
```
