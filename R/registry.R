## ----------------------------------------------------------------------------
## The lens contract & enrolment registry (Stage-0 SS3)
##
## The registry is the durable, versioned record the selector reads. The Stage-1
## battery spans three diversity axes: four inference-family witnesses (dram,
## nimble_gibbs, pg_gibbs, dream_zs), a heavy-tail structure variant (robust_t),
## and a shrinkage spectrum (ridge / lasso / horseshoe) plus bridge_evidence
## diagnostic. Each lens carries its independence key, maturity tier, and (for
## koine-authored lenses) a validation reference; resolve_battery enrols
## family-conditionally and excludes any lens *with a recorded reason* (I1: the
## absence is logged, never hidden).
## ----------------------------------------------------------------------------

#' Construct a lens specification (a registry record)
#'
#' @param id Stable key, e.g. `"dram"`.
#' @param label Human-readable label.
#' @param tier 1 (R/C++ native), 2 (opt-in C++), or 3 (opt-in Python).
#' @param independence_key `list(family =, kernel =)` (Stage-0 SS10).
#' @param role `"opinion"` (enters stacking) or `"diagnostic"` (never stacked).
#' @param exactness `"exact"`, `"approximate"`, or `"resampled"`.
#' @param maturity_tier One of the maturity tiers (Stage-0 SS11).
#' @param family_support Predicate `(family_spec) -> logical`.
#' @param feature_support A `list(random_effects, distributional, nonlinear)`.
#' @param emits_loglik `TRUE` if the lens emits a pointwise log-likelihood
#'   natively; `FALSE` means koine computes it via the family spec.
#' @param fit_fn The adapter `function(pass, lens_spec) -> koine_fb_fit`, or
#'   `NULL` if not yet implemented in this build.
#' @param implemented `TRUE` if `fit_fn` is a working adapter.
#' @param emits_evidence `TRUE` if the lens yields a marginal likelihood.
#' @param structure_variant `TRUE` if the lens fits an alternative model
#'   structure (e.g. a robust-t likelihood) rather than a different inference
#'   method on the same model --- the source of genuine stacking diversity.
#' @param validation_ref Optional identifier (e.g. a ledger cairn id) recording
#'   the SBC / reference-posterior evidence behind a
#'   `koine-implemented-validated` maturity tier (Stage-0 SS11). `NULL` means
#'   no validation evidence is on
#'   record --- a koine-first lens without it is treated as `experimental` by
#'   the maturity gate, never silently trusted (open-arms correctness, SS11).
#' @return A `koine_lens_spec`.
#' @examples
#' lens_spec(
#'   id = "demo", label = "Demonstration lens", tier = 1L,
#'   independence_key = list(family = "random_walk", kernel = "demo"),
#'   role = "opinion", exactness = "exact",
#'   maturity_tier = "reference-validated",
#'   family_support = function(fs) fs$name == "gaussian",
#'   feature_support = list(random_effects = FALSE, distributional = FALSE,
#'                          nonlinear = FALSE),
#'   emits_loglik = FALSE)
#' @export
lens_spec <- function(id, label, tier, independence_key, role, exactness,
                      maturity_tier, family_support, feature_support,
                      emits_loglik, fit_fn = NULL, implemented = FALSE,
                      emits_evidence = FALSE, structure_variant = FALSE,
                      validation_ref = NULL) {
  structure(
    list(id = id, label = label, tier = as.integer(tier),
         independence_key = independence_key, role = role,
         exactness = exactness, maturity_tier = maturity_tier,
         family_support = family_support, feature_support = feature_support,
         emits_loglik = emits_loglik, emits_evidence = emits_evidence,
         fit_fn = fit_fn, implemented = implemented,
         structure_variant = structure_variant,
         validation_ref = validation_ref),
    class = "koine_lens_spec"
  )
}

# Feature support for a plain population-level GLM lens (no RE / dist / nonlin).
.feat_glm_only <- list(random_effects = FALSE, distributional = FALSE,
                       nonlinear = FALSE)

# Feature support for a mixed-model lens (single random intercept).
.feat_re_intercept <- list(random_effects = TRUE, distributional = FALSE,
                           nonlinear = FALSE)

# A heavy-tailed random-effect lens spec at a given Student-t df. The enrolled
# battery is a grid over `nu` (a tail-weight spectrum, like the shrinkage
# spectrum); each `nu` is a distinct structure variant, so stacking can favour
# the random-effect tail-weight the data support. Single source of truth, used
# by both the registry grid and the `fit_glmm_robust_re()` default.
.robust_re_lens_spec <- function(nu) {
  lens_spec(
    id = sprintf("robust_re_t%g", nu),
    label = sprintf("Heavy-tailed RE LMM (Student-t_%g intercepts)", nu),
    tier = 1L,
    independence_key = list(
      family = "random_walk", kernel = sprintf("glmm-robust-re-t%g", nu)),
    role = "opinion", exactness = "approximate",
    maturity_tier = "koine-implemented-validated",
    family_support = function(fs) fs$name == "gaussian",
    feature_support = .feat_re_intercept, emits_loglik = TRUE,
    fit_fn = function(pass, ls) fit_glmm_robust_re(pass, ls, nu = nu),
    implemented = TRUE, structure_variant = TRUE,
    validation_ref = "r-package/2026-06-01-koine-re-stage2b-aghq-robust-re")
}

#' The Stage-1 enrolment registry (leanest cut)
#'
#' @return A list of `koine_lens_spec` records.
#' @examples
#' reg <- koine_stage1_registry()
#' vapply(reg, function(L) L$id, character(1))
#' @export
koine_stage1_registry <- function() {
  list(
    lens_spec(
      id = "dram", label = "Delayed-rejection adaptive Metropolis",
      tier = 1L, independence_key = list(family = "random_walk",
                                         kernel = "dram"),
      role = "opinion", exactness = "exact",
      maturity_tier = "koine-implemented-validated",
      family_support = function(fs) fs$name %in%
        c("gaussian", "binomial", "poisson"),
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_dram, implemented = TRUE,
      validation_ref = "r-package/2026-05-31-koine-lens-sbc-validation"
    ),
    lens_spec(
      id = "nimble_gibbs", label = "NIMBLE Gibbs / Metropolis-within-Gibbs",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "nimble-mwg"),
      role = "opinion", exactness = "exact",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name %in%
        c("gaussian", "binomial", "poisson"),
      feature_support = .feat_glm_only, emits_loglik = TRUE,
      fit_fn = fit_nimble, implemented = TRUE
    ),
    lens_spec(
      id = "pg_gibbs", label = "Polya-Gamma Gibbs (koine-authored, exact)",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "pg-exact"),
      role = "opinion", exactness = "exact",
      maturity_tier = "koine-implemented-validated",
      family_support = function(fs) fs$name == "binomial",
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_pg_gibbs, implemented = TRUE,
      validation_ref = "r-package/2026-05-31-koine-lens-sbc-validation"
    ),
    lens_spec(
      id = "dream_zs", label = "DREAM_ZS differential-evolution MCMC",
      tier = 1L, independence_key = list(family = "population",
                                         kernel = "de-mcmc"),
      role = "opinion", exactness = "approximate",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name %in%
        c("gaussian", "binomial", "poisson"),
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_dream, implemented = TRUE
    ),
    lens_spec(
      id = "robust_t", label = "Robust Student-t linear model (DRAM)",
      tier = 1L, independence_key = list(family = "random_walk",
                                         kernel = "robust-t-dram"),
      role = "opinion", exactness = "exact",
      maturity_tier = "koine-implemented-validated",
      family_support = function(fs) fs$name == "gaussian",
      feature_support = .feat_glm_only, emits_loglik = TRUE,
      fit_fn = fit_robust_t, implemented = TRUE, structure_variant = TRUE,
      validation_ref = "r-package/2026-05-31-koine-lens-sbc-validation"
    ),
    lens_spec(
      id = "horseshoe", label = "Horseshoe sparsity prior (bayesreg Gibbs)",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "bayesreg-hs"),
      role = "opinion", exactness = "exact",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name == "gaussian",
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_horseshoe, implemented = TRUE, structure_variant = TRUE
    ),
    lens_spec(
      id = "lasso", label = "Bayesian lasso (Laplace) prior (bayesreg Gibbs)",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "bayesreg-lasso"),
      role = "opinion", exactness = "exact",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name == "gaussian",
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_lasso, implemented = TRUE, structure_variant = TRUE
    ),
    lens_spec(
      id = "ridge", label = "Bayesian ridge (Gaussian) prior (bayesreg Gibbs)",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "bayesreg-ridge"),
      role = "opinion", exactness = "exact",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name == "gaussian",
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = fit_ridge, implemented = TRUE, structure_variant = TRUE
    ),
    lens_spec(
      id = "glmm_marginal",
      label = "Marginal mixed model (Hamiltonian NUTS, RE integrated out)",
      tier = 1L, independence_key = list(family = "hamiltonian",
                                         kernel = "glmm-marginal-quadrature"),
      role = "opinion", exactness = "approximate",
      maturity_tier = "koine-implemented-validated",
      # Validated (SBC + concordance) for gaussian + binomial, which the
      # Hamiltonian kernel samples cleanly (well-mixed, near-zero divergences).
      # Poisson is NOT enrolled: the NUTS kernel fixes the count-INTERCEPT mixing
      # the random walk could not (the original SBC failure), but a minority of
      # prior-plausible datasets with a large random-intercept SD (tau > 2) give
      # the marginal a variance-component FUNNEL that a static-metric NUTS cannot
      # traverse cheaply -- ~5-10% divergences and order-minutes fits, neither
      # production-viable nor calibration-safe. Marginalising the random effects
      # rules out the usual non-centred fix, so the recorded next step is a
      # funnel-robust metric (Riemannian / Fisher-SoftAbs, reusing the AGHQ
      # information) or a bounded enrolled-tau range. Gaussian + binomial are
      # unaffected.
      family_support = function(fs) fs$name %in% c("gaussian", "binomial"),
      feature_support = .feat_re_intercept, emits_loglik = FALSE,
      fit_fn = fit_glmm_marginal, implemented = TRUE,
      validation_ref = "r-package/2026-05-31-koine-re-lens-sbc-validation"
    ),
    lens_spec(
      id = "nimble_glmm",
      label = "NIMBLE GLMM Gibbs (random effects sampled)",
      tier = 1L, independence_key = list(family = "gibbs",
                                         kernel = "nimble-mwg-glmm"),
      role = "opinion", exactness = "exact",
      maturity_tier = "reference-validated",
      # Gaussian + binomial --- matching the enrolled marginal lens (Poisson is
      # held; see the glmm_marginal record).
      family_support = function(fs) fs$name %in% c("gaussian", "binomial"),
      feature_support = .feat_re_intercept, emits_loglik = FALSE,
      fit_fn = fit_nimble_glmm, implemented = TRUE
    ),
    .robust_re_lens_spec(3),
    .robust_re_lens_spec(8),
    lens_spec(
      id = "bridge_evidence",
      label = "Marginal-likelihood evidence axis (bridge sampling)",
      tier = 1L, independence_key = list(family = "evidence",
                                         kernel = "bridge-sampling"),
      role = "diagnostic", exactness = "approximate",
      maturity_tier = "reference-validated",
      family_support = function(fs) fs$name %in%
        c("gaussian", "binomial", "poisson"),
      feature_support = .feat_glm_only, emits_loglik = FALSE,
      fit_fn = NULL, implemented = TRUE, emits_evidence = TRUE
    )
  )
}

# Fetch a lens record from the default Stage-1 registry by id (used by the lens
# fit functions for their default lens_spec).
.registry_lens <- function(id) {
  reg <- koine_stage1_registry()
  reg[[which(vapply(reg, function(z) z$id, "") == id)]]
}

#' Resolve the battery for a pass (family-conditional enrolment, Stage-0 SS3.3)
#'
#' Enrols a lens iff: its tier is allowed; its family predicate passes; its
#' feature support covers the formula; the log-likelihood is producible (native
#' or via the family spec); the adapter is implemented in this build; and, for a
#' koine-authored lens, its validation evidence is on record (the open-arms
#' correctness gate, SS11 --- an unvalidated first-mover abstains). A lens that
#' fails any test is excluded *with a recorded reason* --- an honest, logged
#' absence, not an error (I1).
#'
#' @param pass A `koine_pass`.
#' @param registry A registry list; defaults to [koine_stage1_registry()].
#' @return A list with `enrolled` (lens specs) and `excluded`
#'   (`list(id =, reason =)` records).
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(30), x = rnorm(30))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' resolved <- resolve_battery(pass)
#' vapply(resolved$enrolled, function(L) L$id, character(1))
#' @export
resolve_battery <- function(pass, registry = koine_stage1_registry()) {
  allowed <- pass$resources$allowed_tiers
  needs_ll <- "log_likelihood" %in% pass$required_outputs
  enrolled <- list()
  excluded <- list()
  for (L in registry) {
    # Diagnostic-role records (e.g. the evidence axis) are not opinion lenses:
    # they are computed post-hoc from the opinion fits, never enrolled into
    # stacking and never an exclusion error (SS6.1, SS9). Skip them here.
    if (identical(L$role, "diagnostic")) {
      next
    }
    reason <- NULL
    if (!L$tier %in% allowed) {
      reason <- sprintf("tier %d not in allowed_tiers", L$tier)
    } else if (!isTRUE(L$family_support(pass$family))) {
      reason <- sprintf("family '%s' unsupported by lens", pass$family$name)
    } else if (isTRUE(pass$has_re) &&
                 !isTRUE(L$feature_support$random_effects)) {
      # A population-level lens cannot honour a random-effects formula; it is
      # excluded with a recorded reason rather than silently fitting the wrong
      # (fixed-effects-only) model (I1).
      reason <- "feature_unsupported: lens has no random-effects support"
    } else if (!isTRUE(pass$has_re) &&
                 isTRUE(L$feature_support$random_effects)) {
      # A mixed-model lens needs a random-effects term; on a population pass it
      # is the population lenses that apply.
      reason <- "feature_unsupported: lens requires a random-effects term"
    } else if (isTRUE(pass$has_re) &&
                 isTRUE(L$feature_support$random_effects) &&
                 !isTRUE(pass$re_supported)) {
      # The random-effects structure (random slopes, multiple or crossed/nested
      # groupings) is outside this build's single-random-intercept support.
      reason <- paste0("feature_unsupported: only a single random intercept ",
                       "on one grouping factor is supported")
    } else if (needs_ll && !L$emits_loglik &&
                 is.null(pass$family$loglik_pointwise)) {
      reason <- "log_likelihood required but not producible"
    } else if (!isTRUE(L$implemented)) {
      reason <- "lens not yet implemented in this build"
    } else if (.lens_needs_validation(L) &&
                 (is.null(L$validation_ref) ||
                    !nzchar(as.character(L$validation_ref)))) {
      # Open-arms correctness gate (SS11): a koine-authored lens may not
      # contribute to triangulation until its validation evidence is on
      # record. An unvalidated first-mover abstains, never silently trusted.
      reason <- "validation_gate: koine-first lens not yet validated"
    }
    if (is.null(reason)) {
      enrolled[[length(enrolled) + 1L]] <- L
    } else {
      excluded[[length(excluded) + 1L]] <- list(id = L$id, reason = reason)
    }
  }
  list(enrolled = enrolled, excluded = excluded)
}

#' @export
print.koine_lens_spec <- function(x, ...) {
  cat("<koine lens_spec '", x$id, "'>\n", sep = "")
  cat("  ", x$label, "\n", sep = "")
  cat("  tier ", x$tier, " | role ", x$role, " | ",
      x$independence_key$family, "/", x$independence_key$kernel,
      " | ", x$maturity_tier, "\n", sep = "")
  cat("  implemented: ", x$implemented, "\n", sep = "")
  invisible(x)
}
