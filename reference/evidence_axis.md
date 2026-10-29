# Evidence-axis diagnostic (Stage-0 SS9.3)

Per-lens log marginal likelihood (model evidence) via bridge sampling
(Meng & Wong 1996; Gronau et al. 2017), as a Bayes-factor cross-check on
the battery. This axis has `role = "diagnostic"`: it is a scalar per
lens and **never enters stacking** (the synthesis assembly excludes
diagnostic-role output). The estimate is **prior-conditional** — it is
computed under koine's canonical reference prior (coefficients N(0,
`beta_prior_sd`); Gaussian dispersion via a log-sigma N(0,
`log_sigma_prior_sd`) prior), and Bayes factors are sensitive to that
prior, so the result is flagged accordingly. Lenses that fit an
alternative model structure (a different likelihood, e.g. robust-t) are
excluded from the cross-check: their evidence is not comparable on the
shared prior.

## Usage

``` r
evidence_axis(
  fb_fits,
  pass,
  beta_prior_sd = 5,
  log_sigma_prior_sd = 2,
  log_tau_prior_sd = 2
)
```

## Arguments

- fb_fits:

  A list of contributed `koine_fb_fit` objects.

- pass:

  The originating `koine_pass`.

- beta_prior_sd:

  Prior SD on the coefficients (default 5, matching the lenses'
  weakly-informative default).

- log_sigma_prior_sd:

  Prior SD on log-sigma for Gaussian families (default 2).

- log_tau_prior_sd:

  Prior SD on log-tau for mixed-model passes (default 2, matching the RE
  lenses' weakly-informative default).

## Value

A list with `available` (logical), `prior_conditional` (always TRUE),
and, when available, `evidence` (a data frame of `lens_id`, `log_ml`,
`log_ml_se`) plus `note`. `role` is documented as `"diagnostic"`.

## Details

Requires the suggested bridgesampling package; absent it, returns a
structured unavailable record rather than failing (I1).

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40))
pass <- koine_pass(y ~ x, d, gaussian())
fit <- fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
evidence_axis(list(fit), pass)$evidence
#>   lens_id    log_ml    log_ml_se
#> 1    dram -60.36968 0.0001461722
# }
```
