# Package index

## Request and contract

The engine-neutral request, the evaluable family adapter, and the
uniform lens-output contract.

- [`koine_pass()`](https://max578.github.io/koine/reference/koine_pass.md)
  :

  Construct a koine `pass`

- [`koine_can_contribute()`](https://max578.github.io/koine/reference/koine_can_contribute.md)
  : Can koine contribute a diverse opinion for this model? (pre-flight)

- [`as_family_spec()`](https://max578.github.io/koine/reference/as_family_spec.md)
  :

  Adapt a base-R family to an evaluable `family_spec`

- [`family_spec()`](https://max578.github.io/koine/reference/family_spec.md)
  : Construct an evaluable family specification

- [`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md) :

  Construct an `fb_fit`

- [`validate_fb_fit()`](https://max578.github.io/koine/reference/validate_fb_fit.md)
  :

  Validate an `fb_fit` against the field rules

- [`fb_fit_diagnostics()`](https://max578.github.io/koine/reference/fb_fit_diagnostics.md)
  :

  Recover convergence diagnostics from an `fb_fit`'s draws

- [`provenance_record()`](https://max578.github.io/koine/reference/provenance_record.md)
  : Build a provenance record

- [`validate_provenance()`](https://max578.github.io/koine/reference/validate_provenance.md)
  : Validate a provenance record's mandatory tuple

## Battery and enrolment registry

- [`lens_spec()`](https://max578.github.io/koine/reference/lens_spec.md)
  : Construct a lens specification (a registry record)
- [`koine_stage1_registry()`](https://max578.github.io/koine/reference/koine_stage1_registry.md)
  : The Stage-1 enrolment registry (leanest cut)
- [`resolve_battery()`](https://max578.github.io/koine/reference/resolve_battery.md)
  : Resolve the battery for a pass (family-conditional enrolment,
  Stage-0 SS3.3)

## Opinion lenses

Inference-family witnesses and model-structure variants.

- [`fit_dram()`](https://max578.github.io/koine/reference/fit_dram.md) :
  Fit the DRAM lens
- [`fit_nimble()`](https://max578.github.io/koine/reference/fit_nimble.md)
  : Fit the NIMBLE-Gibbs lens
- [`fit_pg_gibbs()`](https://max578.github.io/koine/reference/fit_pg_gibbs.md)
  : Fit the Polya-Gamma Gibbs lens (logistic regression)
- [`fit_dream()`](https://max578.github.io/koine/reference/fit_dream.md)
  : Fit the DREAM_ZS lens
- [`fit_robust_t()`](https://max578.github.io/koine/reference/fit_robust_t.md)
  : Fit the robust-t lens (Student-t linear model)
- [`fit_ridge()`](https://max578.github.io/koine/reference/fit_ridge.md)
  : Fit the Bayesian ridge lens
- [`fit_lasso()`](https://max578.github.io/koine/reference/fit_lasso.md)
  : Fit the Bayesian lasso lens
- [`fit_horseshoe()`](https://max578.github.io/koine/reference/fit_horseshoe.md)
  : Fit the horseshoe (sparsity) lens

## Mixed-model lenses

Random-effects lenses — a single random intercept and the
breeder-standard multi-environment-trial (MET) structure of two
crossed/nested intercepts. The marginal lens integrates the random
effects out; the NIMBLE lens samples them; the robust-RE lens gives them
heavy tails.

- [`fit_glmm_marginal()`](https://max578.github.io/koine/reference/fit_glmm_marginal.md)
  : Fit the marginal mixed-model lens

- [`fit_nimble_glmm()`](https://max578.github.io/koine/reference/fit_nimble_glmm.md)
  : Fit the NIMBLE GLMM lens

- [`fit_glmm_robust_re()`](https://max578.github.io/koine/reference/fit_glmm_robust_re.md)
  : Fit the heavy-tailed random-effect lens (Student-t random
  intercepts)

- [`koine_pass_from_design()`](https://max578.github.io/koine/reference/koine_pass_from_design.md)
  :

  Construct a koine `pass` from design arrays (the IR-targetable
  builder)

## Synthesis and gate

- [`synthesise()`](https://max578.github.io/koine/reference/synthesise.md)
  : Synthesise the fourth opinion for a pass (Stage-0 SS1 flow)
- [`koine_opinion_draws()`](https://max578.github.io/koine/reference/koine_opinion_draws.md)
  : koine's opinion as a posterior draws array (the stacking mixture)
- [`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)
  : Assemble per-lens pointwise log-likelihood matrices (Stage-0 SS6.1)
- [`gate_psis()`](https://max578.github.io/koine/reference/gate_psis.md)
  : The two-tier gate (Stage-0 SS7) — tier 1 (PSIS reliability)
- [`stack_opinion()`](https://max578.github.io/koine/reference/stack_opinion.md)
  : LOO-stacking of opinion lenses (Stage-0 SS6.2)

## Diagnostic layer

- [`weight_concentration()`](https://max578.github.io/koine/reference/weight_concentration.md)
  : Weight-concentration diagnostic (Stage-0 SS9)
- [`linear_pool_envelope()`](https://max578.github.io/koine/reference/linear_pool_envelope.md)
  : Linear-pool envelope diagnostic (Stage-0 SS9.2)
- [`evidence_axis()`](https://max578.github.io/koine/reference/evidence_axis.md)
  : Evidence-axis diagnostic (Stage-0 SS9.3)
- [`independence_report()`](https://max578.github.io/koine/reference/independence_report.md)
  : Build the independence report for the contributing lenses (Stage-0
  SS10)

## Validation and maturity gate

- [`sbc_check()`](https://max578.github.io/koine/reference/sbc_check.md)
  : Simulation-based calibration of a lens (Stage-0 SS11; Talts et al.
  2018)
- [`reference_posterior_check()`](https://max578.github.io/koine/reference/reference_posterior_check.md)
  : Exact grid-based reference-posterior check (Stage-0 SS11)
- [`gate_corroboration()`](https://max578.github.io/koine/reference/gate_corroboration.md)
  : The corroboration / maturity gate (Stage-0 SS11)
- [`lens_effective_maturity()`](https://max578.github.io/koine/reference/lens_effective_maturity.md)
  : Effective maturity tier of a lens after the validation gate
- [`koine_lens_validation_setup()`](https://max578.github.io/koine/reference/koine_lens_validation_setup.md)
  : Build the validation setup for a koine-authored lens

## Abstention and interchange

- [`abstention()`](https://max578.github.io/koine/reference/abstention.md)
  : Construct a structured abstention

- [`is_abstention()`](https://max578.github.io/koine/reference/is_abstention.md)
  : Is this object an abstention?

- [`read_idata_nc()`](https://max578.github.io/koine/reference/read_idata_nc.md)
  :

  Read an ArviZ DataTree NetCDF-4 file into an `fb_fit`
