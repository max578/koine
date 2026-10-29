# Construct a lens specification (a registry record)

Construct a lens specification (a registry record)

## Usage

``` r
lens_spec(
  id,
  label,
  tier,
  independence_key,
  role,
  exactness,
  maturity_tier,
  family_support,
  feature_support,
  emits_loglik,
  fit_fn = NULL,
  implemented = FALSE,
  emits_evidence = FALSE,
  structure_variant = FALSE,
  validation_ref = NULL
)
```

## Arguments

- id:

  Stable key, e.g. `"dram"`.

- label:

  Human-readable label.

- tier:

  1 (R/C++ native), 2 (opt-in C++), or 3 (opt-in Python).

- independence_key:

  `list(family =, kernel =)` (Stage-0 SS10).

- role:

  `"opinion"` (enters stacking) or `"diagnostic"` (never stacked).

- exactness:

  `"exact"`, `"approximate"`, or `"resampled"`.

- maturity_tier:

  One of the maturity tiers (Stage-0 SS11).

- family_support:

  Predicate `(family_spec) -> logical`.

- feature_support:

  A `list(random_effects, distributional, nonlinear, met)`. `met`
  (default `FALSE`) flags a mixed-model lens that also covers the MET
  structure of two crossed/nested random intercepts.

- emits_loglik:

  `TRUE` if the lens emits a pointwise log-likelihood natively; `FALSE`
  means koine computes it via the family spec.

- fit_fn:

  The adapter `function(pass, lens_spec) -> koine_fb_fit`, or `NULL` if
  not yet implemented in this build.

- implemented:

  `TRUE` if `fit_fn` is a working adapter.

- emits_evidence:

  `TRUE` if the lens yields a marginal likelihood.

- structure_variant:

  `TRUE` if the lens fits an alternative model structure (e.g. a
  robust-t likelihood) rather than a different inference method on the
  same model — the source of genuine stacking diversity.

- validation_ref:

  Optional identifier (e.g. a ledger cairn id) recording the SBC /
  reference-posterior evidence behind a `koine-implemented-validated`
  maturity tier (Stage-0 SS11). `NULL` means no validation evidence is
  on record — a koine-first lens without it is treated as `experimental`
  by the maturity gate, never silently trusted (open-arms correctness,
  SS11).

## Value

A `koine_lens_spec`.

## Examples

``` r
lens_spec(
  id = "demo", label = "Demonstration lens", tier = 1L,
  independence_key = list(family = "random_walk", kernel = "demo"),
  role = "opinion", exactness = "exact",
  maturity_tier = "reference-validated",
  family_support = function(fs) fs$name == "gaussian",
  feature_support = list(random_effects = FALSE, distributional = FALSE,
                         nonlinear = FALSE, met = FALSE),
  emits_loglik = FALSE)
#> <koine lens_spec 'demo'>
#>   Demonstration lens
#>   tier 1 | role opinion | random_walk/demo | reference-validated
#>   implemented: FALSE
```
