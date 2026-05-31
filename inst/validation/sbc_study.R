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

# Auto-run only when executed as a script (Rscript), not when source()d --- when
# sourced, sys.nframe() > 0 and the caller invokes run_validation_study() itself.
if (sys.nframe() == 0L) {
  print(run_validation_study())
}
