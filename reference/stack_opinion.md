# LOO-stacking of opinion lenses (Stage-0 SS6.2)

Computes stacking weights over the retained lenses and the stacked
leave-one-out predictive density (the stacking objective evaluated at
the optimum). Weights are retained for the diagnostic surface; they are
never shipped as a pooled posterior (I2).

## Usage

``` r
stack_opinion(loo_list)
```

## Arguments

- loo_list:

  A named list of loo objects (the gate's `loo`).

## Value

A list with `weights` (named simplex), `elpd_stacked`, and
`n_components`.

## References

Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking
to average Bayesian predictive distributions. *Bayesian Analysis*,
13(3), 917-1007.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = 1 + rnorm(60), x = rnorm(60))
pass <- koine_pass(y ~ x, d, gaussian())
fits <- list(
  fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L),
  fit_robust_t(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L))
gate <- gate_psis(assemble_loglik(fits, pass))
if (length(gate$kept) >= 2L) stack_opinion(gate$loo)
#> $weights
#>      dram  robust_t 
#> 0.6685778 0.3314222 
#> 
#> $elpd_stacked
#> [1] -78.43638
#> 
#> $n_components
#> [1] 2
#> 
# }
```
