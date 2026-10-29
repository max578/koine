# Validate a provenance record's mandatory tuple

Validate a provenance record's mandatory tuple

## Usage

``` r
validate_provenance(x)
```

## Arguments

- x:

  A `koine_provenance` object.

## Value

`x`, invisibly, if every mandatory field is present and non-`NULL`;
otherwise raises `error.provenance_incomplete`.

## Examples

``` r
prov <- provenance_record(
  method = "dram", method_version = "0.1.0", impl_language = "R",
  seed = 1L, derived_seed = 42L, ir_hash = "abc123",
  maturity_tier = "reference-validated",
  independence_key = list(family = "random_walk", kernel = "dram"),
  role = "opinion", exactness = "exact")
validate_provenance(prov)
```
