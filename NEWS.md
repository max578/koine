# koine 0.0.0.9000

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
