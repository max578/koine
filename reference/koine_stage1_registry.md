# The Stage-1 enrolment registry (leanest cut)

The Stage-1 enrolment registry (leanest cut)

## Usage

``` r
koine_stage1_registry()
```

## Value

A list of `koine_lens_spec` records.

## Examples

``` r
reg <- koine_stage1_registry()
vapply(reg, function(L) L$id, character(1))
#>  [1] "dram"            "nimble_gibbs"    "pg_gibbs"        "dream_zs"       
#>  [5] "robust_t"        "horseshoe"       "lasso"           "ridge"          
#>  [9] "glmm_marginal"   "nimble_glmm"     "robust_re_t3"    "robust_re_t8"   
#> [13] "bridge_evidence"
```
