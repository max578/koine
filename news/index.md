# Changelog

## koine (development version)

- The package has a logo.
- README lists related packages and carries check, coverage and
  r-universe badges.
- CI runs lintr (report-only) and test coverage.

## koine 0.3.1

- The package vignette, now titled *A second opinion for a mixed model*,
  was restructured to the orchestra’s vignette quality bar: a stated
  question, a fixed Why/What/Do/Read/Limits/What-to-read-next/Reproduce
  shape, every fitted example shrunk to render in under two minutes,
  captioned colourblind-safe `ggplot2` figures in place of base
  graphics, tables through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), and a
  stated limit that the random-effects lenses’ fixed `N(0, 5)`
  coefficient prior dominates on raw-scale (tens- or hundreds-scale)
  responses unless the caller standardises the response and covariates
  before building the pass.
- **Fix**:
  [`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)
  no longer silently drops a lens that abstains at fit time. Previously
  the fit-time `koine_abstention` a lens returns had no `$status` field,
  so the “contributed” filter dropped it without recording why, and the
  synthesis-level abstention then blamed the PSIS gate
  (`psis_unreliable`) even when no lens ever reached it. The lens’s real
  reason (e.g. `convergence_failed`, `feature_unsupported`) now
  propagates to the synthesis-level `reason` when it is the sole
  distinct reason among the declining lenses, and every declining lens
  is carried into `lenses_dropped` with that reason.
- `inst/INTEGRATION.md` refreshed to describe the MET random-effects
  structure (two crossed/nested intercepts) as in-corner, add the
  `groups` argument to the IR-mapping table, and correct the stale
  version/licence status block.
- Corrected a licence-metadata contradiction: `README.md` and
  `CITATION.cff` previously stated no licence was declared and the
  package was undistributed, while `DESCRIPTION` has declared
  `MIT + file LICENSE` since 2026-06-08 and the package is pushed to
  `max578/koine`. Both now state MIT and the correct
  installation/contributing instructions.
- Corrected the RE-unsupported synthesis abstention’s `detail` text,
  which described the supported corner as “a single random intercept on
  one grouping factor” – stale since the MET (two crossed/nested
  intercepts) structure landed at 0.3.0. It now names both supported
  shapes.
- Corrected a stale sentence in the *koine: the synthesised fourth
  opinion* vignette that described the Poisson GLMM lens as “held
  experimental pending a better-mixing sampler”; it has been enrolled
  (with a post-fit convergence gate) since 0.2.0.
- Typeset the MET marginal-covariance equation in the vignette as LaTeX
  maths (previously inline code) and added two captioned, interpreted
  figures: the two-lens stacking-weight simplex with `n_eff` annotated,
  and the over-dispersed linear-pool envelope against the fitted data.
- Replaced the em dash (`---`) with the house en dash (`--`) throughout
  `README.md` and the vignette prose.
- Added `@references` (Vehtari, Simpson, Gelman, Yao & Gabry 2024 for
  the PSIS gate; Yao, Vehtari, Simpson & Gelman 2018 for the stacking
  objective) and an `@examples` block for
  [`is_abstention()`](https://max578.github.io/koine/reference/is_abstention.md).
- Fixed the new metadata-consistency test (README/`CITATION.cff` licence
  wording) failing with an error, rather than skipping, under
  `R CMD check`, where the source-tree files it reads are not installed
  alongside the package.
- Corrected `CONTRIBUTING.md`’s stale instruction to branch off a `dev`
  branch that does not exist in this repository; it now says `main`.
- Added `inst/CITATION`, so `citation("koine")` returns a reference.

## koine 0.3.0

### New features

- **The breeder-standard multi-environment-trial (MET) random-effects
  structure is now supported**: two crossed or nested random
  *intercepts*, a variety effect and a trial (environment) effect,
  `(1 | variety) + (1 | trial)`. koine now corroborates a production MET
  model *directly* — on the variety-plus-trial structure the production
  model actually fits — rather than on a single-random-effect reduction
  that discards the trial structure. Whether the two groupings are
  crossed or nested is resolved from the data, never from the formula.
- The MET marginal stays **exact for the Gaussian LMM**: with two
  crossed groupings the marginal covariance is dense, so the
  single-random-intercept per-group quadrature no longer applies, but
  the marginal is still the exact multivariate normal
  `y ~ N(X beta, sigma^2 I + tau_variety^2 Z_variety Z_variety' + tau_trial^2 Z_trial Z_trial')`.
  [`fit_glmm_marginal()`](https://max578.github.io/koine/reference/fit_glmm_marginal.md)
  samples `(beta, sigma, tau_variety, tau_trial)` on that exact marginal
  (a dense Cholesky of `V`) with a coordinate-wise slice kernel,
  integrating the random effects out where the sealed engines sample or
  Laplace-approximate them. The exact marginal log-likelihood was
  validated to machine precision against an independent `lme4` ML fit
  (the Independent Oracle Principle), and the MET fit concords with the
  NIMBLE MET reference (which *samples* the two crossed random effects).
  The draws now carry one `tau_<grouping>` per random factor.
- The lift is **additive**. The single-random-intercept path `(1 | g)`
  is unchanged (byte-identical code path and accessors); the MET
  structure is the new capability, dispatched only when two groupings
  are present. No exported function gained or changed a required
  argument, so the public API and the C1/C3 triangulation contract hold.
  [`koine_pass_from_design()`](https://max578.github.io/koine/reference/koine_pass_from_design.md)
  gained a `groups` argument (a named list of two grouping vectors)
  alongside the existing `group`.
- A non-Gaussian MET — whose dense crossed-GLMM marginal has no exact
  quadrature — returns a lens-scope abstention (`feature_unsupported`)
  rather than an ungrounded approximation; the single-random-intercept
  GLMM path (Gaussian, binomial, Poisson) is untouched. The heavy-tailed
  single-RE lens
  ([`fit_glmm_robust_re()`](https://max578.github.io/koine/reference/fit_glmm_robust_re.md))
  is excluded from a MET pass with a recorded reason rather than
  dropping the second grouping.

### Minor improvements and fixes

- The closed abstention taxonomy gained two fact-grounding reasons,
  `data_ungrounded` and `provenance_unverified` (Independent Oracle
  Principle). A converged, high-ESS lens fit on an ungrounded input
  still carries that input’s ungroundedness; these let a lens decline on
  a *fact* failure (an input marked ungrounded in the orchestra
  provenance vocabulary), not only on the existing numerical failures
  (convergence, calibration, ESS).
- **Fixed a latent NIMBLE compile failure for an intercept-only fixed
  part** (`K = 1`). NIMBLE’s `inprod()` over a length-1 vector generates
  non-compiling Eigen C++, so the NIMBLE GLMM lenses now write the
  `K = 1` fixed-effects mean as the scalar `beta[1]` and keep `inprod()`
  for `K >= 2`. The bug never fired before because every test used a
  `K = 2` fixed part; the MET model (a commonly intercept-only fixed
  part) surfaced it.

## koine 0.2.0

- **The Poisson GLMM (single random intercept, log link) is now
  ENROLLED** for the marginal mixed-model lens (`fit_glmm_marginal`) and
  the NIMBLE reference lens (`fit_nimble_glmm`), validated by
  simulation-based calibration (SBC clean over the contributing
  simulations: every parameter’s rank uniform, p \> 0.16) and
  cross-engine concordance against NIMBLE (median standardised
  discrepancy 0.23, SD ratios 0.94–1.17) on the
  identifiable-variance-component regime.
- The Poisson marginal is sampled by a funnel-safe **blocked slice
  kernel**: `log tau` is slice-sampled, and the coefficients are
  slice-sampled along their `tau`-resolved principal axes (the scaled
  eigenvectors of the Laplace covariance of the coefficient
  conditional). This clears the variance-component funnel neck where a
  static-metric Hamiltonian kernel funnels (~50% divergences, `ess_tau`
  ~ 3) and a Laplace independence proposal freezes; the Gaussian and
  binomial marginals keep the joint NUTS kernel. A value-only marginal
  log-posterior (`.glmm_marginal_logpost_value`) backs the slice kernel
  (it skips the score- identity gradient quadrature the slice never
  uses).
- At a large random-intercept SD (`tau` beyond ~2, the ~8% upper tail of
  the prior) the marginalised funnel is too weakly identified for any
  marginal sampler to mix — with the random effects integrated out, the
  non-centred remedy is unavailable. A **post-fit convergence gate**
  makes the lens ABSTAIN (`convergence_failed`, a new first-class
  abstention reason — too few effective draws or chains that disagree)
  rather than contribute a badly-mixed opinion. koine contributes a
  calibrated Poisson opinion where the variance component is
  identifiable, and abstains where it is not (a first-class outcome,
  SS3).

## koine 0.1.1

- [`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md)
  now mixes each parameter only over the lenses for which it is
  comparable: the linear-predictor coefficients over every contributing
  lens, but a scale parameter (`sigma`, `tau`) only over the
  family-likelihood lenses. A likelihood-variant lens (the robust-t
  lens, which carries its own Student-t likelihood) reports its scale on
  a different metric (a Student-t scale rather than a Gaussian SD), so
  it is excluded from the scale mixture — and flagged in the new
  `koine_scale_restricted` attribute — while still contributing to the
  coefficients. This fixes a synthesised-`sigma` bias when a
  likelihood-variant lens carried stacking weight (it dragged `sigma`
  toward the t-scale); the fixed effects are unchanged. The draws are a
  per-parameter marginal representation, which is what a
  parameter-marginal triangulation needs.

## koine 0.1.0

First tagged release — the version carrying the IR-targetable,
re-entrant model builder
([`koine_pass_from_design()`](https://max578.github.io/koine/reference/koine_pass_from_design.md))
and the integration surface a parent triangulation panel drives, so
koine clears the backend contract’s construction, output, and
determinism clauses. Licence declared **MIT**.

### Stage-2 — random-effects support (single random intercept)

- [`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md)
  now parses an `lme4`-style random intercept `(1 | g)` on one grouping
  factor into a fixed-effects design plus a random-effects design;
  richer structure (random slopes, multiple or crossed/nested groupings)
  is parsed and flagged unsupported so the battery abstains
  (`feature_unsupported`) rather than silently fitting the wrong model.
  This also closes a latent gap whereby a random-effects formula
  previously slipped through as a wrongly parsed population model.
- [`resolve_battery()`](https://max578.github.io/koine/reference/resolve_battery.md)
  gained the feature-support gate (a lens enrols only when its
  random-effects support matches the pass), with directional recorded
  reasons.
- Two orthogonal mixed-model lenses, stacked on the identical
  per-observation *integrated* predictive (the random effect
  marginalised under each draw’s variance components, so PSIS stays in
  its reliable regime):
  [`fit_glmm_marginal()`](https://max578.github.io/koine/reference/fit_glmm_marginal.md)
  — a koine-authored No-U-Turn (Hamiltonian) sampler that integrates the
  random effects out (exact for the Gaussian LMM, adaptive Gauss-Hermite
  quadrature for a GLMM) and rides the marginal’s analytic gradient (the
  score identity for a GLMM, closed form for the Gaussian) — and
  [`fit_nimble_glmm()`](https://max578.github.io/koine/reference/fit_nimble_glmm.md)
  — NIMBLE Gibbs that samples the random effects. The two carry
  identical priors, so their agreement is method-independent consensus.
  Enrolled (validated) for the Gaussian LMM and the binomial GLMM (the
  Poisson GLMM is implemented but held experimental — see below).
- A heavy-tailed random-effect lens,
  [`fit_glmm_robust_re()`](https://max578.github.io/koine/reference/fit_glmm_robust_re.md):
  a Gaussian LMM whose random intercepts are Student-t (`b_g ~ t_nu`),
  integrated out via the Normal scale-mixture (Gauss-Laguerre over the
  mixing scale). It is enrolled as a **spectrum over the degrees of
  freedom** (a grid `nu in {3, 8}`, analogous to the shrinkage
  spectrum), so stacking can favour the random-effect tail-weight the
  data support. It is the within-battery analogue of the population
  robust-t lens — where the Gaussian-RE lenses corroborate (weights
  collapse), the t-RE lenses earn real stacking weight when group
  effects are heavy-tailed (`n_eff ~ 1.75` on `t_3` effects), so koine’s
  RE opinion becomes sensitive to random-effect structure rather than a
  weighted echo of one assumption.
- Adaptive Gauss-Hermite quadrature for the non-Gaussian random-effect
  marginal (per-group Newton mode-finding + curvature-scaled nodes,
  closed-form canonical-link scores added to `family_spec`). Fixed
  Gauss-Hermite is non-convergent for the sharply-peaked Poisson count
  integrand; adaptive Gauss-Hermite is accurate to ~1e-3 across the
  prior range and makes the binomial GLMM validate more accurately.
- Family-appropriate weakly-informative priors (`.re_prior()`): the
  log-link Poisson uses a tighter coefficient / random-intercept-SD
  prior than the Gaussian and binomial families, whose `N(0, 5)` prior
  is fine.
- [`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)
  and
  [`evidence_axis()`](https://max578.github.io/koine/reference/evidence_axis.md)
  gained random-effects branches built on a single shared
  marginalisation core.
- A NaN-safe sampler guard: `.dram_one_chain()` and the marginal
  log-posteriors guard both Metropolis acceptance ratios against
  non-finite log-posteriors, so an extreme proposal is rejected (and a
  chain started from an impossible point escapes) rather than halting on
  a `NaN` acceptance test. Finite cases are unchanged.
- The validation harness covers every mixed-model lens; the enrolled
  ones are validated by SBC plus cross-implementation concordance
  against NIMBLE (recorded study in `inst/validation/`). The Poisson
  GLMM remains held experimental and is not enrolled. The Hamiltonian
  kernel fixes the count-intercept mixing that the random walk could not
  (the original failure), but a minority of prior-plausible datasets
  with a large random-intercept SD give the marginal a
  variance-component funnel that a static-metric sampler cannot traverse
  cheaply (occasional divergences and order-minutes fits). Because the
  random effects are marginalised, the usual non-centred remedy does not
  apply; the recorded next step is a funnel-robust metric (Riemannian /
  Fisher information, reused from the quadrature) or a bounded
  enrolled-tau range.
- [`koine_pass_from_design()`](https://max578.github.io/koine/reference/koine_pass_from_design.md)
  — an IR-targetable, re-entrant, deparse-free model builder that
  constructs a pass directly from a named design matrix, a response, a
  family, and an optional grouping vector (no model formula, no
  [`match.call()`](https://rdrr.io/r/base/match.call.html) capture, no
  session-global state). This is the programmatic entry a parent
  triangulation panel drives (the backend contract’s construction
  clause); it preserves non-syntactic canonical column names such as
  `(Intercept)` verbatim.
- Integration surface for a parent triangulation panel (the contract is
  in `inst/INTEGRATION.md`):
  [`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md)
  materialises the single synthesised opinion as a posterior
  `draws_array` over the parameters its contributing lenses share,
  sampling each lens in proportion to its stacking weight (a mixture
  over koine’s own battery, with the stacking weights and stacked
  predictive elpd surfaced as attributes — never a pool across the
  parent’s independent engines).
  [`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)
  gains `keep_fits` to retain the per-lens fits for it, and
  [`koine_can_contribute()`](https://max578.github.io/koine/reference/koine_can_contribute.md)
  is a cheap pre-flight capability predicate for the parent’s
  backend-availability gate.

### Stage-1 — the synthesis spine

First development scaffold — the Stage-1 increment-1 spine.

- The engine-neutral request
  ([`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md))
  with content-hashed provenance, the evaluable family adapter
  ([`as_family_spec()`](https://max578.github.io/koine/reference/as_family_spec.md),
  Gaussian / Binomial / Poisson), and the uniform lens-output contract
  ([`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md))
  with a complete, validated provenance record.

- The enrolment registry and family-conditional
  [`resolve_battery()`](https://max578.github.io/koine/reference/resolve_battery.md)
  for the leanest Stage-1 cut (DRAM, NIMBLE-Gibbs, Polya-Gamma Gibbs,
  DREAM); only the DRAM lens carries a working adapter in this
  increment, the other three are registered but excluded with a recorded
  reason.

- The first opinion lens: a self-contained, koine-authored
  delayed-rejection adaptive Metropolis sampler
  ([`fit_dram()`](https://max578.github.io/koine/reference/fit_dram.md)).

- The synthesis spine: pointwise log-likelihood assembly
  ([`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)),
  the tier-1 PSIS reliability gate
  ([`gate_psis()`](https://max578.github.io/koine/reference/gate_psis.md)),
  predictive-space LOO-stacking
  ([`stack_opinion()`](https://max578.github.io/koine/reference/stack_opinion.md)),
  and the top-level
  [`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)
  — which returns a first-class abstention when fewer than two
  independent opinion lenses survive.

- The dual-reader NetCDF interchange
  ([`read_idata_nc()`](https://max578.github.io/koine/reference/read_idata_nc.md),
  RNetCDF canonical with an ncdf4 fallback), enforcing the four hard
  read assertions.

- Abstention
  ([`abstention()`](https://max578.github.io/koine/reference/abstention.md))
  and the typed error taxonomy as first-class outcomes.

- The full Stage-1 battery across three diversity axes: four
  inference-family witnesses (DRAM, NIMBLE-Gibbs, Polya-Gamma Gibbs,
  DREAM_ZS), a heavy-tail structure variant
  ([`fit_robust_t()`](https://max578.github.io/koine/reference/fit_robust_t.md)),
  and a shrinkage spectrum
  ([`fit_ridge()`](https://max578.github.io/koine/reference/fit_ridge.md)
  /
  [`fit_lasso()`](https://max578.github.io/koine/reference/fit_lasso.md)
  /
  [`fit_horseshoe()`](https://max578.github.io/koine/reference/fit_horseshoe.md)),
  with the weight-concentration diagnostic and the `weight_degenerate`
  abstention.

- The diagnostic layer is complete: alongside weight concentration,
  [`linear_pool_envelope()`](https://max578.github.io/koine/reference/linear_pool_envelope.md)
  reports an over-dispersed reference bracket and flags disagreeing
  lenses leave-one-out, and
  [`evidence_axis()`](https://max578.github.io/koine/reference/evidence_axis.md)
  reports per-lens marginal likelihoods (bridge sampling,
  prior-conditional). Both carry `role = "diagnostic"` and are surfaced
  by
  [`synthesise()`](https://max578.github.io/koine/reference/synthesise.md),
  never pooled into the opinion.

- A validation / maturity gate for koine-authored lenses (the open-arms
  correctness obligation): simulation-based calibration
  ([`sbc_check()`](https://max578.github.io/koine/reference/sbc_check.md)),
  an exact grid-based reference-posterior check
  ([`reference_posterior_check()`](https://max578.github.io/koine/reference/reference_posterior_check.md)),
  and the corroboration gate
  ([`gate_corroboration()`](https://max578.github.io/koine/reference/gate_corroboration.md)).
  [`resolve_battery()`](https://max578.github.io/koine/reference/resolve_battery.md)
  now excludes a koine-first lens until its validation evidence is on
  record — an unvalidated first-mover abstains. The DRAM, Polya-Gamma
  Gibbs, and robust-t lenses are validated (the recorded study is in
  `inst/validation/`).

- Release-grade documentation and tooling: a runnable `@examples` block
  on every exported function, an `en-GB` spell-check gate
  (`inst/WORDLIST`), a `pkgdown` site, a `README`, a `CITATION.cff`, an
  `air.toml` / `.lintr` style configuration, and an R-CMD-check GitHub
  Actions workflow. `R CMD check --as-cran` is clean (no ERRORs or
  WARNINGs).

No licence is declared yet (pending the required consultation); the
package is local-only and not distributed.
