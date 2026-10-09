# koine lens validation (Stage-0 SS11)

koine's *open-arms* posture lets it implement a method before its authors ship a
reference implementation. A koine-authored lens therefore vouches for its own
correctness, which weakens the "agreement across independent implementations is
evidence" story. To keep the fourth opinion sound, every koine-authored lens
must clear an assessment gate before it may contribute to triangulation, and it
carries a maturity tier:

```
reference-validated  >  koine-implemented-validated  >  experimental
```

A lens with no validation evidence on record is treated as `experimental` by the
corroboration gate (`gate_corroboration()`) and **abstains** — it is never
silently trusted.

## What this study runs

`sbc_study.R` runs, for each koine-authored lens (`dram`, `pg_gibbs`,
`robust_t`):

1. **Simulation-based calibration** (SBC; Talts et al. 2018). True parameters are
   drawn from the lens's prior, a dataset is simulated from the lens's *own*
   generative model (a Student-t likelihood for `robust_t`, a Bernoulli-logit
   model for `pg_gibbs`, a Gaussian model for `dram`), the lens is fitted, and
   the rank of each true parameter among the thinned posterior draws is recorded.
   A correctly calibrated sampler yields ranks uniform on `0:n_keep`; a
   per-parameter chi-squared test (Bonferroni-adjusted across parameters) detects
   departures.

2. **Reference-posterior agreement** (`reference_posterior_check()`). On a sample
   of simulated datasets the lens posterior is compared to an exact reference
   obtained by dense grid integration of the same unnormalised log-posterior.
   The grid reference is independent of any sampler and needs no external engine,
   so it is runnable in this standalone, pre-graduation package. Discrepancy is
   measured in reference standard deviations.

A lens is judged `koine-implemented-validated` when SBC passes **and** the median
reference-posterior discrepancy stays within tolerance.

## Scope boundary — cross-engine concordance (deferred)

SS11 names a third criterion, cross-engine concordance (closeness to the
independent-engine consensus). It requires the engine panel that koine joins at
graduation and is **not** run here. This is an explicit, recorded scope boundary,
not a silent omission: a koine-authored lens is `koine-implemented-validated` on
the SBC + reference-posterior evidence and is **discounted** by the independence
registry accordingly (`gate_corroboration()` returns a 0.5 weight for
koine-validated lenses). The third criterion is scheduled for the graduation
milestone.

## Results

The recorded verdicts live in `results/`:

- `results/summary.csv` / `results/summary.rds` — one row per lens (SBC minimum
  p-value, SBC pass/fail, reference-posterior median discrepancy, verdict).
- `results/sbc_<lens>.rds` — the full per-lens record (rank matrix, uniformity
  table, reference discrepancies, verdict).

### Random-effects lenses (`run_re_validation_study()`)

The mixed-model lenses (`glmm_marginal`, `nimble_glmm`, the `robust_re_t<nu>`
spectrum) are validated by **SBC plus cross-implementation concordance against
NIMBLE** (which samples the random effects; koine's marginal lens integrates them
out), recorded in `results/summary_re.{csv,rds}` and
`results/sbc_glmm_marginal_<family>.rds`. The `glmm_marginal` lens is validated
for **gaussian, binomial, and poisson**.

The **Poisson** marginal rides a funnel-safe blocked double-slice kernel (slice
`log tau`; slice the coefficients along their `tau`-resolved Laplace principal
axes). At a large random-intercept SD (`tau` beyond ~2, the upper tail of the
prior) the marginalised variance-component funnel is too weakly identified for any
marginal sampler to mix, so a post-fit convergence gate makes the lens **abstain**
(`convergence_failed`) rather than contribute a badly-mixed opinion. The recorded
Poisson study is therefore over the *contributing* sims (71/80 here; 8 abstained,
1 capped): SBC clean (every parameter's rank uniform, min p = 0.16), concordance
median standardised discrepancy 0.23 vs NIMBLE. Because the Poisson gate fits the
tail fully before abstaining, the recorded result is generated at lighter settings
than the default (`fit_control = list(n_iter = 1200L, n_burn = 600L, n_chains =
2L)`, `n_sims = 80L`, `n_keep = 50L`, `n_conc = 3L`).

## Reproduce

```r
# after install
Rscript inst/validation/sbc_study.R
# in development
devtools::load_all(); source("inst/validation/sbc_study.R"); run_validation_study()
```

The study is seeded; re-running reproduces the recorded ranks and verdicts.
