# Synthesise the fourth opinion for a pass (Stage-0 SS1 flow)

Resolves the battery, fits each opinion lens, assembles pointwise
log-likelihoods, applies the PSIS gate, stacks the survivors, and
reports a weight-concentration diagnostic — or returns a first-class
[`abstention()`](https://max578.github.io/koine/reference/abstention.md)
when fewer than two independent opinion lenses survive. The result is
either a `koine_result` or a `koine_abstention`, never a silent partial.

## Usage

``` r
synthesise(
  pass,
  registry = koine_stage1_registry(),
  on_degenerate = c("flag", "abstain"),
  keep_fits = FALSE
)
```

## Arguments

- pass:

  A `koine_pass`.

- registry:

  A registry list; defaults to
  [`koine_stage1_registry()`](https://max578.github.io/koine/reference/koine_stage1_registry.md).

- on_degenerate:

  What to do when the stacking weights are degenerate (the battery
  corroborated rather than diversified): `"flag"` (default) returns the
  opinion with the degeneracy recorded in its diagnostics; `"abstain"`
  returns a `weight_degenerate` abstention (SS8).

- keep_fits:

  If `TRUE`, attach the contributing per-lens fits to the result (its
  `fits` field), so
  [`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md)
  can materialise the stacking-mixture posterior without re-fitting.
  Default `FALSE` (the diagnostic-only contract is unchanged).

## Value

A `koine_result` or a `koine_abstention`.

## Examples

``` r
d <- data.frame(y = rnorm(30), x = rnorm(30))
pass <- koine_pass(y ~ x, d, gaussian())
# with a single-lens battery, koine abstains rather than guessing
res <- synthesise(pass, registry = list(koine_stage1_registry()[[1]]))
is_abstention(res)
#> [1] TRUE
```
