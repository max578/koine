# Build the validation setup for a koine-authored lens

Returns the prior draw, the model simulator, the sampler call, the
parameters to rank, and a grid-exact reference-posterior builder for one
of koine's koine-authored lenses (`"dram"`, `"pg_gibbs"`, `"robust_t"`).
This is the shared harness behind
[`sbc_check()`](https://max578.github.io/koine/reference/sbc_check.md)
and
[`reference_posterior_check()`](https://max578.github.io/koine/reference/reference_posterior_check.md)
for the recorded validation study; SBC simulates from the model returned
here, so the calibration test sees the lens's own
`{prior, model, sampler}` triple.

## Usage

``` r
koine_lens_validation_setup(
  lens_id,
  n_obs = 40L,
  n_pred = 1L,
  nu = 4,
  fit_control = list(n_iter = 1500L, n_burn = 750L, n_chains = 2L),
  family = c("gaussian", "binomial", "poisson"),
  n_groups = 12L,
  n_per_group = 6L
)
```

## Arguments

- lens_id:

  One of `"dram"`, `"pg_gibbs"`, `"robust_t"` (population) or
  `"glmm_marginal"`, `"nimble_glmm"` (mixed-model lenses).

- n_obs:

  Observations per simulated dataset (population lenses).

- n_pred:

  Number of predictors beyond the intercept.

- nu:

  Student-t degrees of freedom for `"robust_t"` (the generative model
  and the sampler share this value).

- fit_control:

  A list of sampler controls passed to the lens fit.

- family:

  For the mixed-model lenses, the response family (`"gaussian"`,
  `"binomial"`, or `"poisson"`); ignored by the population lenses, which
  carry their own family.

- n_groups:

  Number of random-effect groups (mixed-model lenses).

- n_per_group:

  Observations per group (mixed-model lenses).

## Value

A list with `generate`, `fit_fn`, `param_names`, and `reference`.

## Examples

``` r
setup <- koine_lens_validation_setup("dram")
setup$param_names
#> [1] "b_(Intercept)" "b_x1"          "sigma"        
```
