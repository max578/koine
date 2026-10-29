# Build a provenance record

Build a provenance record

## Usage

``` r
provenance_record(
  method,
  method_version,
  impl_language,
  seed,
  derived_seed,
  ir_hash,
  maturity_tier,
  independence_key,
  role,
  exactness,
  timestamp_iso8601 = NULL,
  lib_versions = NULL,
  device = NULL,
  py_pins = NULL,
  trace = NULL,
  structure_variant = FALSE
)
```

## Arguments

- method:

  Lens method id (e.g. `"dram"`).

- method_version:

  Version string of the method implementation.

- impl_language:

  `"R"`, `"C++"`, or `"Python"`.

- seed:

  Master seed from the pass.

- derived_seed:

  Per-lens derived seed.

- ir_hash:

  Content hash of the originating pass.

- maturity_tier:

  One of the maturity tiers (Stage-0 SS11).

- independence_key:

  `list(family =, kernel =)` (Stage-0 SS10).

- role:

  `"opinion"` or `"diagnostic"`.

- exactness:

  `"exact"`, `"approximate"`, or `"resampled"`.

- timestamp_iso8601:

  ISO-8601 timestamp; defaults to now (UTC).

- lib_versions:

  Optional named list of library versions.

- device:

  Optional device string.

- py_pins:

  Optional named list of Python pins (Tier-3 only).

- trace:

  Optional list of per-stage outcomes.

- structure_variant:

  `TRUE` if the lens fits an alternative model structure (a different
  likelihood or prior) rather than the canonical reference posterior.
  Carried so the diagnostic layer (the shared-prior evidence axis,
  SS9.3) and the independence registry's `specification` axis (SS10.1)
  can tell structure variants from reference-posterior lenses.

## Value

A `koine_provenance` object.

## Examples

``` r
provenance_record(
  method = "dram", method_version = "0.1.0", impl_language = "R",
  seed = 1L, derived_seed = 42L, ir_hash = "abc123",
  maturity_tier = "reference-validated",
  independence_key = list(family = "random_walk", kernel = "dram"),
  role = "opinion", exactness = "exact")
```
