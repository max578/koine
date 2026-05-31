# koine 0.0.0.9000

## Stage-2 --- random-effects support (single random intercept)

* `koine_pass()` now parses an `lme4`-style random intercept `(1 | g)` on one
  grouping factor into a fixed-effects design plus a random-effects design;
  richer structure (random slopes, multiple or crossed/nested groupings) is
  parsed and flagged unsupported so the battery abstains
  (`feature_unsupported`) rather than silently fitting the wrong model. This
  also closes a latent gap whereby a random-effects formula previously slipped
  through as a wrongly parsed population model.
* `resolve_battery()` gained the feature-support gate (a lens enrols only when
  its random-effects support matches the pass), with directional recorded
  reasons.
* Two orthogonal mixed-model lenses, stacked on the identical per-observation
  *integrated* predictive (the random effect marginalised under each draw's
  variance components, so PSIS stays in its reliable regime): `fit_glmm_marginal()`
  --- a koine-authored random-walk sampler that integrates the random effects
  out (exact for the Gaussian LMM, Gauss-Hermite quadrature for a GLMM) --- and
  `fit_nimble_glmm()` --- NIMBLE Gibbs that samples the random effects. The two
  carry identical priors, so their agreement is method-independent consensus.
  Enrolled (validated) for the Gaussian LMM and the binomial GLMM. The Poisson
  GLMM quadrature is implemented and tested for correctness but held
  experimental: under the default coefficient prior its prior-predictive is
  pathological, so calibration (SBC) cannot yet be established --- it is not
  enrolled (see the validation study).
* `assemble_loglik()` and `evidence_axis()` gained random-effects branches built
  on a single shared marginalisation core (the per-observation integrated
  predictive and the per-group marginal likelihood).
* The validation harness covers the mixed-model lenses (`koine_lens_validation_setup()`
  family / grouping arguments); the marginal lens is validated by SBC plus
  cross-implementation concordance against NIMBLE for the Gaussian and binomial
  families (recorded study in `inst/validation/`).

## Stage-1 --- the synthesis spine

First development scaffold --- the Stage-1 increment-1 spine.

* The engine-neutral request (`koine_pass()`) with content-hashed provenance,
  the evaluable family adapter (`as_family_spec()`, Gaussian / Binomial /
  Poisson), and the uniform lens-output contract (`fb_fit()`) with a complete,
  validated provenance record.
* The enrolment registry and family-conditional `resolve_battery()` for the
  leanest Stage-1 cut (DRAM, NIMBLE-Gibbs, Polya-Gamma Gibbs, DREAM); only the
  DRAM lens carries a working adapter in this increment, the other three are
  registered but excluded with a recorded reason.
* The first opinion lens: a self-contained, koine-authored delayed-rejection
  adaptive Metropolis sampler (`fit_dram()`).
* The synthesis spine: pointwise log-likelihood assembly (`assemble_loglik()`),
  the tier-1 PSIS reliability gate (`gate_psis()`), predictive-space
  LOO-stacking (`stack_opinion()`), and the top-level `synthesise()` --- which
  returns a first-class abstention when fewer than two independent opinion
  lenses survive.
* The dual-reader NetCDF interchange (`read_idata_nc()`, RNetCDF canonical with
  an ncdf4 fallback), enforcing the four hard read assertions.
* Abstention (`abstention()`) and the typed error taxonomy as first-class,
  honest outcomes.
* The full Stage-1 battery across three diversity axes: four inference-family
  witnesses (DRAM, NIMBLE-Gibbs, Polya-Gamma Gibbs, DREAM_ZS), a heavy-tail
  structure variant (`fit_robust_t()`), and a shrinkage spectrum
  (`fit_ridge()` / `fit_lasso()` / `fit_horseshoe()`), with the
  weight-concentration diagnostic and the `weight_degenerate` abstention.
* The diagnostic layer is complete: alongside weight concentration,
  `linear_pool_envelope()` reports an over-dispersed honest reference bracket
  and flags disagreeing lenses leave-one-out, and `evidence_axis()` reports
  per-lens marginal likelihoods (bridge sampling, prior-conditional). Both carry
  `role = "diagnostic"` and are surfaced by `synthesise()`, never pooled into
  the opinion.
* A validation / maturity gate for koine-authored lenses (the open-arms
  correctness obligation): simulation-based calibration (`sbc_check()`), an
  exact grid-based reference-posterior check (`reference_posterior_check()`), and
  the corroboration gate (`gate_corroboration()`). `resolve_battery()` now
  excludes a koine-first lens until its validation evidence is on record --- an
  unvalidated first-mover abstains. The DRAM, Polya-Gamma Gibbs, and robust-t
  lenses are validated (the recorded study is in `inst/validation/`).

* Release-grade documentation and tooling: a runnable `@examples` block on every
  exported function, an `en-GB` spell-check gate (`inst/WORDLIST`), a `pkgdown`
  site, a `README`, a `CITATION.cff`, an `air.toml` / `.lintr` style
  configuration, and an R-CMD-check GitHub Actions workflow. `R CMD check
  --as-cran` is clean (no ERRORs or WARNINGs).

No licence is declared yet (pending the required consultation); the package is
local-only and not distributed.
