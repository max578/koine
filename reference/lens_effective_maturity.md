# Effective maturity tier of a lens after the validation gate

Resolves a lens spec's declared maturity against its `validation_ref`: a
koine-first lens without recorded validation evidence is demoted to
`experimental` (it may not corroborate). Reference-validated lenses and
validation-backed koine lenses keep their declared tier.

## Usage

``` r
lens_effective_maturity(lens_spec)
```

## Arguments

- lens_spec:

  A `koine_lens_spec`.

## Value

The effective maturity tier (a string).

## Examples

``` r
lens_effective_maturity(koine_stage1_registry()[[1]])
#> [1] "koine-implemented-validated"
```
