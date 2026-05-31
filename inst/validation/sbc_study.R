## ----------------------------------------------------------------------------
## Recorded validation study for koine's koine-authored lenses (Stage-0 SS11).
##
## For each koine-first lens (dram, pg_gibbs, robust_t) this runs:
##   * simulation-based calibration (SBC; Talts et al. 2018) against the lens's
##     own generative model, and
##   * a grid-exact reference-posterior check on a sample of simulated datasets.
##
## A lens is judged `koine-implemented-validated` when SBC ranks are uniform
## (every parameter's Bonferroni-adjusted chi-squared p-value clears alpha) AND
## the reference-posterior discrepancy stays within tolerance on the sample.
## Cross-engine concordance --- the third SS11 criterion --- needs the engine
## panel and is an explicit scope boundary deferred to graduation (NOT run here).
##
## Reproduce after install with:  Rscript inst/validation/sbc_study.R
## In development:                devtools::load_all(); source(this file)
## ----------------------------------------------------------------------------

if (!"koine" %in% loadedNamespaces()) {
  if (requireNamespace("koine", quietly = TRUE)) {
    library(koine)
  } else {
    stop("install koine or run devtools::load_all() before this script.")
  }
}

run_validation_study <- function(out_dir = NULL, n_sims = 160L, n_keep = 64L,
                                  n_ref = 16L, alpha = 0.05, verbose = TRUE) {
  if (is.null(out_dir)) {
    out_dir <- system.file("validation", package = "koine")
    if (!nzchar(out_dir)) out_dir <- file.path("inst", "validation")
  }
  res_dir <- file.path(out_dir, "results")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  lenses <- list(
    dram = list(n_obs = 50L, n_pred = 1L,
                ctrl = list(n_iter = 4000L, n_burn = 2000L, n_chains = 3L)),
    pg_gibbs = list(n_obs = 200L, n_pred = 1L,
                    ctrl = list(n_iter = 3000L, n_burn = 1500L, n_chains = 3L)),
    robust_t = list(n_obs = 50L, n_pred = 1L,
                    ctrl = list(n_iter = 4000L, n_burn = 2000L, n_chains = 3L))
  )

  summary_rows <- list()
  for (id in names(lenses)) {
    cfg <- lenses[[id]]
    if (verbose) message("== validating ", id, " ==")
    setup <- koine_lens_validation_setup(
      id, n_obs = cfg$n_obs, n_pred = cfg$n_pred, fit_control = cfg$ctrl)

    sbc <- sbc_check(setup$generate, setup$fit_fn, setup$param_names,
                     n_sims = n_sims, n_keep = n_keep, n_bins = 8L,
                     alpha = alpha, seed = 20260531L)

    # reference-posterior discrepancy over a sample of simulated datasets.
    disc <- vapply(seq_len(n_ref), function(i) {
      g <- setup$generate(900000L + i)
      fit <- setup$fit_fn(g$pass)
      if (is_abstention(fit)) return(NA_real_)
      ref <- setup$reference(g)
      chk <- reference_posterior_check(fit, ref$log_post_unnorm, ref$grid,
                                       tol = 0.2)
      if (!isTRUE(chk$applicable)) return(NA_real_)
      max(chk$discrepancy)
    }, numeric(1))
    disc <- disc[is.finite(disc)]
    ref_pass <- length(disc) > 0L && stats::median(disc) < 0.2

    verdict <- if (sbc$passed && ref_pass) {
      "koine-implemented-validated"
    } else {
      "experimental"
    }
    out <- list(lens = id, config = cfg, sbc = sbc,
                reference_discrepancy = disc, ref_pass = ref_pass,
                verdict = verdict)
    saveRDS(out, file.path(res_dir, paste0("sbc_", id, ".rds")))
    summary_rows[[id]] <- data.frame(
      lens = id, n_eff_sims = sbc$n_eff_sims, sbc_min_p = min(sbc$uniformity$p_value),
      sbc_passed = sbc$passed, ref_median_disc = stats::median(disc),
      ref_passed = ref_pass, verdict = verdict, stringsAsFactors = FALSE)
    if (verbose) {
      message(sprintf("   SBC min p = %.3f (pass=%s); ref median disc = %.3f (pass=%s) -> %s",
                      min(sbc$uniformity$p_value), sbc$passed,
                      stats::median(disc), ref_pass, verdict))
    }
  }

  summ <- do.call(rbind, summary_rows)
  saveRDS(summ, file.path(res_dir, "summary.rds"))
  utils::write.csv(summ, file.path(res_dir, "summary.csv"), row.names = FALSE)
  summ
}

## ----------------------------------------------------------------------------
## Mixed-model (random-intercept) lens validation (Stage-2).
##
## glmm_marginal is koine-authored, so it must clear the gate before it may
## contribute. Two criteria are used:
##   * Simulation-based calibration (SBC; Talts et al. 2018) on (beta, sigma,
##     tau) --- AUTHORITATIVE, including the variance component. Data are
##     simulated from the full generative model (draw tau, b_g ~ N(0, tau^2),
##     y), and the lens targets the marginal posterior; rank uniformity on
##     (beta, sigma, tau) is the correct criterion because marginalising b is
##     exact.
##   * Cross-implementation concordance against nimble_glmm --- NIMBLE samples
##     the random effects (Gibbs) and is reference-validated; the marginal lens
##     integrates them out (quadrature). Two independent implementations of the
##     same marginal posterior must agree to within posterior noise. This stands
##     in for the grid-exact reference-posterior check used for the population
##     lenses: that check (a uniform-grid Riemann integral) is reliable for a
##     low-dimensional population posterior but numerically fragile for the
##     marginal mixed-model posterior (the weakly-identified variance component
##     makes the grid reference's SD collapse intermittently), so it is NOT used
##     as an RE criterion --- recorded here explicitly, never silently dropped.
## Cross-ENGINE concordance (the three sealed engines) remains the third SS11
## criterion, deferred to graduation as for every koine lens.
## ----------------------------------------------------------------------------

run_re_validation_study <- function(out_dir = NULL, n_sims = 120L,
                                    n_keep = 64L, n_conc = 6L, alpha = 0.05,
                                    conc_tol = 0.5,
                                    families = c("gaussian", "binomial",
                                                 "poisson"),
                                    verbose = TRUE) {
  if (is.null(out_dir)) {
    out_dir <- system.file("validation", package = "koine")
    if (!nzchar(out_dir)) out_dir <- file.path("inst", "validation")
  }
  res_dir <- file.path(out_dir, "results")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  ctrl <- list(n_iter = 4000L, n_burn = 2000L, n_chains = 3L)
  have_nimble <- requireNamespace("nimble", quietly = TRUE)

  # Cross-implementation concordance: do the marginal (quadrature) and NIMBLE
  # (sampled-b) lenses agree on the marginal posterior? Standardised by NIMBLE's
  # posterior SD; the max over (beta, sigma, tau), median over datasets.
  concordance <- function(fam) {
    if (!have_nimble) {
      return(list(available = FALSE, disc = NA_real_, pass = NA))
    }
    sm <- koine_lens_validation_setup("glmm_marginal", family = fam,
                                      n_groups = 15L, n_per_group = 6L,
                                      fit_control = ctrl)
    disc <- vapply(seq_len(n_conc), function(i) {
      g <- sm$generate(800000L + i)
      fm <- sm$fit_fn(g$pass)
      fn <- fit_nimble_glmm(g$pass, n_iter = ctrl$n_iter, n_burn = ctrl$n_burn,
                            n_chains = ctrl$n_chains)
      if (is_abstention(fm) || is_abstention(fn)) return(NA_real_)
      a <- posterior::summarise_draws(fm$draws, mean = mean)
      b <- posterior::summarise_draws(fn$draws, mean = mean, sd = stats::sd)
      max(abs(a$mean - b$mean) / pmax(b$sd, 1e-6))
    }, numeric(1))
    disc <- disc[is.finite(disc)]
    list(available = TRUE, disc = stats::median(disc),
         pass = length(disc) > 0L && stats::median(disc) < conc_tol)
  }

  rows <- list()
  for (fam in families) {
    if (verbose) message("== glmm_marginal / ", fam, " ==")
    setup <- koine_lens_validation_setup("glmm_marginal", family = fam,
                                         n_groups = 15L, n_per_group = 6L,
                                         fit_control = ctrl)
    sbc <- sbc_check(setup$generate, setup$fit_fn, setup$param_names,
                     n_sims = n_sims, n_keep = n_keep, n_bins = 8L,
                     alpha = alpha, seed = 20260531L)
    conc <- concordance(fam)
    # Gate: SBC is authoritative; concordance must also hold where checkable.
    verdict <- if (sbc$passed &&
                     (!isTRUE(conc$available) || isTRUE(conc$pass))) {
      "koine-implemented-validated"
    } else {
      "experimental"
    }
    out <- list(lens = "glmm_marginal", family = fam, sbc = sbc,
                concordance = conc, verdict = verdict)
    saveRDS(out, file.path(res_dir, paste0("sbc_glmm_marginal_", fam, ".rds")))
    rows[[length(rows) + 1L]] <- data.frame(
      lens = "glmm_marginal", family = fam, n_eff_sims = sbc$n_eff_sims,
      sbc_min_p = min(sbc$uniformity$p_value), sbc_passed = sbc$passed,
      concordance_median = conc$disc, concordance_pass = conc$pass,
      verdict = verdict, stringsAsFactors = FALSE)
    if (verbose) {
      message(sprintf(
        "   SBC min p = %.3f (pass=%s); concordance median = %.3f (pass=%s) -> %s",
        min(sbc$uniformity$p_value), sbc$passed, conc$disc, conc$pass, verdict))
    }
  }
  summ <- do.call(rbind, rows)
  saveRDS(summ, file.path(res_dir, "summary_re.rds"))
  utils::write.csv(summ, file.path(res_dir, "summary_re.csv"), row.names = FALSE)
  summ
}

# Auto-run only when executed as a script (Rscript), not when source()d --- when
# sourced, sys.nframe() > 0 and the caller invokes run_validation_study() itself.
if (sys.nframe() == 0L) {
  print(run_validation_study())
  print(run_re_validation_study())
}
