# Weight-concentration diagnostic (Stage-0 SS9)

Summarises how concentrated the stacking weights are. Concentration on
one component means the battery added no predictive diversity (the
lenses corroborate rather than diversify); spread weights mean genuine
model-structure diversity was realised. Reports the effective number of
components `n_eff = 1 / sum(w^2)`, the maximum weight, the entropy, and
a `degenerate` flag.

## Usage

``` r
weight_concentration(weights, n_eff_threshold = 1.5)
```

## Arguments

- weights:

  A named numeric vector of stacking weights.

- n_eff_threshold:

  Below this effective component count the weights are flagged
  degenerate. Default 1.5.

## Value

A list with `weights`, `max_weight`, `n_eff`, `entropy`, and
`degenerate`.

## Examples

``` r
# spread weights: genuine diversity; concentrated weights: corroboration
weight_concentration(c(dram = 0.55, robust_t = 0.45))
#> $weights
#>     dram robust_t 
#>     0.55     0.45 
#> 
#> $max_weight
#> [1] 0.55
#> 
#> $n_eff
#> [1] 1.980198
#> 
#> $entropy
#> [1] 0.6881388
#> 
#> $degenerate
#> [1] FALSE
#> 
weight_concentration(c(dram = 0.97, robust_t = 0.03))
#> $weights
#>     dram robust_t 
#>     0.97     0.03 
#> 
#> $max_weight
#> [1] 0.97
#> 
#> $n_eff
#> [1] 1.061797
#> 
#> $entropy
#> [1] 0.1347422
#> 
#> $degenerate
#> [1] TRUE
#> 
```
