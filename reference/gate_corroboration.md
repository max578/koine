# The corroboration / maturity gate (Stage-0 SS11)

Decides whether a lens of a given maturity tier may corroborate the
fourth opinion, and at what weight discount. A koine-first lens
(`koine-implemented-validated`) may corroborate only when validation
evidence is on record (`validation_ref`); without it the lens is treated
as `experimental` and may not corroborate — an unvalidated first-mover
is never silently trusted (open-arms correctness, SS6/SS11).

## Usage

``` r
gate_corroboration(maturity_tier, validation_ref = NULL)
```

## Arguments

- maturity_tier:

  One of `"reference-validated"`, `"koine-implemented-validated"`,
  `"experimental"`.

- validation_ref:

  Optional identifier of the validation evidence (e.g. a ledger cairn
  id). Required to back the `koine-implemented-validated` tier.

## Value

A list with `effective_tier`, `may_corroborate` (logical), and
`discount` (1 reference-validated, 0.5 koine-validated, 0 experimental).

## Examples

``` r
gate_corroboration("reference-validated")
#> $effective_tier
#> [1] "reference-validated"
#> 
#> $may_corroborate
#> [1] TRUE
#> 
#> $discount
#> [1] 1
#> 
gate_corroboration("koine-implemented-validated", validation_ref = "cairn-1")
#> $effective_tier
#> [1] "koine-implemented-validated"
#> 
#> $may_corroborate
#> [1] TRUE
#> 
#> $discount
#> [1] 0.5
#> 
gate_corroboration("koine-implemented-validated") # unbacked -> experimental
#> $effective_tier
#> [1] "experimental"
#> 
#> $may_corroborate
#> [1] FALSE
#> 
#> $discount
#> [1] 0
#> 
```
