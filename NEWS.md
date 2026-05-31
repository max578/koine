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
  Enrolled (validated) for the Gaussian LMM and the binomial GLMM.
* A third, heavy-tailed random-effect lens, `fit_glmm_robust_re()`: a Gaussian
  LMM whose random intercepts are Student-t (`b_g ~ t_nu`), integrated out via
  the Normal scale-mixture (Gauss-Laguerre over the mixing scale). It is the
  within-battery analogue of the population robust-t lens --- where the two
  Gaussian-RE lenses corroborate (weights collapse), the t-RE lens earns real
  stacking weight when group effects are heavy-tailed, so koine's RE opinion
  becomes sensitive to random-effect structure rather than a weighted echo of
  one assumption.
* Adaptive Gauss-Hermite quadrature for the non-Gaussian random-effect marginal
  (per-group Newton mode-finding + curvature-scaled nodes, closed-form
  canonical-link scores added to `family_spec`). Fixed Gauss-Hermite is
  non-convergent for the sharply-peaked Poisson count integrand; adaptive
  Gauss-Hermite is accurate to ~1e-3 across the prior range and makes the
  binomial GLMM validate more accurately.
* Family-appropriate weakly-informative priors (`.re_prior()`): the log-link
  Poisson uses a tighter coefficient / random-intercept-SD prior than the
  Gaussian and binomial families, whose `N(0, 5)` prior is fine.
* `assemble_loglik()` and `evidence_axis()` gained random-effects branches built
  on a single shared marginalisation core.
* The validation harness covers all three mixed-model lenses; they are validated
  by SBC plus cross-implementation concordance against NIMBLE (recorded study in
  `inst/validation/`). The Poisson GLMM remains held experimental and is not
  enrolled: adaptive Gauss-Hermite fixed its marginal-likelihood quadrature (its
  variance
  component and slope now calibrate), but the random-walk sampler mixes the
  strongly-confounded count intercept too poorly for SBC to pass --- a
  better-mixing sampler is the recorded next step.

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
