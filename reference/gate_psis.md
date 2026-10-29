# The two-tier gate (Stage-0 SS7) — tier 1 (PSIS reliability)

Computes per-lens PSIS-LOO with a sample-size-dependent Pareto-k
threshold `min(1 - 1 / log10(S), 0.7)`. Lenses whose worst k exceeds the
threshold are dropped with a recorded reason (never silently). This
build implements the drop action; K-fold refit is a later stage.

## Usage

``` r
gate_psis(loglik_list)
```

## Arguments

- loglik_list:

  Output of
  [`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md).

## Value

A list with `loo` (named list of loo objects for retained lenses),
`kept` (ids), `dropped` (`list(id =, reason =, k =)`), and `threshold`.

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72), 1-58.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40))
pass <- koine_pass(y ~ x, d, gaussian())
fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
gate_psis(assemble_loglik(list(fit), pass))
#> $loo
#> list()
#> 
#> $kept
#> NULL
#> 
#> $dropped
#> $dropped[[1]]
#> $dropped[[1]]$id
#> [1] "dram"
#> 
#> $dropped[[1]]$reason
#> [1] "psis_unreliable"
#> 
#> $dropped[[1]]$k
#> [1] 0.7647691
#> 
#> 
#> 
#> $threshold
#> [1] 0.6400484
#> 
# }
```
