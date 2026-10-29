# Build the independence report for the contributing lenses (Stage-0 SS10)

One row per contributing opinion lens, carrying its two-level key,
maturity tier, and a corroboration discount. The discount reflects:
koine-first maturity (SS11), within-family kernel redundancy, and
membership of a collapse grouping. This feeds the downstream
fourth-vs-three grading.

## Usage

``` r
independence_report(fb_fits)
```

## Arguments

- fb_fits:

  A list of contributed `koine_fb_fit` objects.

## Value

A data frame with columns `lens_id`, `family`, `kernel`,
`maturity_tier`, `discount`, and `collapsed_witness`.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40))
pass <- koine_pass(y ~ x, d, gaussian())
fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
independence_report(list(fit))
#>   lens_id      family kernel               maturity_tier discount
#> 1    dram random_walk   dram koine-implemented-validated      0.5
#>   collapsed_witness
#> 1              <NA>
# }
```
