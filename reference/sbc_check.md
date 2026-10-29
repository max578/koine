# Simulation-based calibration of a lens (Stage-0 SS11; Talts et al. 2018)

Repeatedly draws true parameters from the prior, simulates a dataset
from the lens's own generative model, fits the lens, and records the
rank of each true parameter among the (thinned) posterior draws. If the
lens samples its posterior correctly, the ranks are uniform on
`0:n_keep`; systematic departures (a chi-squared test per parameter)
reveal a miscalibrated sampler.

## Usage

``` r
sbc_check(
  generate,
  fit_fn,
  param_names,
  n_sims = 128L,
  n_keep = 64L,
  n_bins = 8L,
  alpha = 0.05,
  seed = 1L
)
```

## Arguments

- generate:

  A function `(seed) -> list(theta = named-numeric, pass = koine_pass)`:
  draws true parameters from the prior and builds a pass whose data is
  simulated from the model at `theta`. `theta`'s names must match the
  posterior draw variables to be ranked.

- fit_fn:

  A function `(pass) -> koine_fb_fit` (the lens, with its sampler
  controls already set).

- param_names:

  Character vector of parameters to rank (a subset of the names of
  `theta` and of the fit's draw variables).

- n_sims:

  Number of SBC replications.

- n_keep:

  Posterior draws retained per fit (the ranks lie in `0:n_keep`); draws
  are thinned to this many to keep the rank statistic near-independent.

- n_bins:

  Number of histogram bins for the chi-squared uniformity test.

- alpha:

  Family-wise significance level (Bonferroni-adjusted across parameters)
  below which calibration is judged to have failed.

- seed:

  Master seed; per-replication seeds are derived from it.

## Value

A list with `ranks` (an `n_sims` x `length(param_names)` matrix),
`n_keep`, `uniformity` (a data frame of `param`, `chisq`, `df`,
`p_value`), `alpha_adj`, and `passed` (logical: every parameter's
p-value exceeds the adjusted alpha).

## Details

SBC validates the `{prior, model, sampler}` triple as a whole, so the
`generate` callback must simulate data from exactly the model the lens
assumes (e.g. a Student-t likelihood for the robust-t lens), using the
exact prior the lens places.

## Examples

``` r
# \donttest{
setup <- koine_lens_validation_setup(
  "dram", n_obs = 30L,
  fit_control = list(n_iter = 800L, n_burn = 400L, n_chains = 2L))
res <- sbc_check(setup$generate, setup$fit_fn, setup$param_names,
                 n_sims = 16L, n_keep = 16L, n_bins = 4L)
res$passed
#> [1] TRUE
# }
```
