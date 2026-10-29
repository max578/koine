# Exact grid-based reference-posterior check (Stage-0 SS11)

Compares a lens's posterior to an exact reference computed by dense grid
integration of the same unnormalised log-posterior — a reference that is
independent of any sampler and needs no external engine, so it is
runnable in this standalone package. Feasible only for low-dimensional
models (the grid is the product over parameters); for larger models it
returns a structured "not applicable" record rather than a misleading
number (I1).

## Usage

``` r
reference_posterior_check(
  fit,
  log_post_unnorm,
  grid,
  tol = 0.15,
  max_cells = 2e+05
)
```

## Arguments

- fit:

  A contributed `koine_fb_fit` from the lens under test.

- log_post_unnorm:

  A function `(theta) -> numeric`: the lens's unnormalised log-posterior
  (likelihood + prior) at a named parameter vector on the **draw scale**
  (the same variable names the fit's draws carry).

- grid:

  A named list of numeric grid vectors, one per parameter, naming the
  parameters in draw-scale order.

- tol:

  Maximum tolerated absolute standardised discrepancy between the lens
  and reference posterior mean (in reference SDs) for a pass.

- max_cells:

  Refuse (return not-applicable) above this many grid cells.

## Value

A list with `applicable` (logical), and when applicable `reference` and
`lens` summaries (per-parameter mean and SD), `discrepancy` (per
parameter, in reference SDs), and `passed`.

## Examples

``` r
# \donttest{
setup <- koine_lens_validation_setup(
  "dram", n_obs = 30L,
  fit_control = list(n_iter = 1500L, n_burn = 750L, n_chains = 2L))
g <- setup$generate(1L)
fit <- setup$fit_fn(g$pass)
ref <- setup$reference(g, n_pts = 21L)
reference_posterior_check(fit, ref$log_post_unnorm, ref$grid)$passed
#> [1] TRUE
# }
```
