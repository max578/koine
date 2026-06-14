## ----------------------------------------------------------------------------
## The diagnostic layer (Stage-0 SS9) --- surfaced, never pooled
##
## Three named, first-class diagnostics, none ever shipped as the opinion:
##   1. engine-weight concentration --- `weight_concentration()` (synthesise.R),
##      the in-stack degeneracy check; applied here to the contributing lenses.
##   2. linear-pool envelope        --- `linear_pool_envelope()`, an
##      over-dispersed honest reference bracket; flags lenses that disagree.
##   3. evidence axis               --- `evidence_axis()`, per-lens marginal
##      likelihood (bridge sampling); role = "diagnostic", NEVER stacked.
##
## The Stage-0 contract phrases (1) and (2) over the three sealed engines. koine
## is a standalone, pre-graduation package with no engines of its own, so the
## diagnostics are built engine-agnostically: they operate over whatever
## predictive components they are handed. In standalone use the caller supplies
## koine's contributing opinion lenses; once koine joins the engine panel the
## same machinery takes the three engines' predictives. This keeps SS9's intent
## (a disagreement signal, an over-dispersed bracket, an evidence cross-check)
## without hard-wiring "three engines".
## ----------------------------------------------------------------------------

# Posterior of the conditional mean mu_i = linkinv(eta_i) per draw, for a fit.
# Returns an S x N matrix on the response (mean) scale, or NULL if the fit's
# draws do not expose the expected b_<col> coefficients (e.g. a diagnostic-role
# fit). The mean scale is family-agnostic and comparable across every lens that
# shares the linear-predictor structure --- so the envelope brackets the
# model-structure disagreement, not the likelihood-tail differences.
.fit_mean_response <- function(fit, design, family) {
  if (identical(fit$status, "abstained") ||
        identical(fit$provenance$role, "diagnostic")) {
    return(NULL)
  }
  mat <- unclass(posterior::as_draws_matrix(fit$draws))
  bcols <- paste0("b_", colnames(design$X))
  if (!all(bcols %in% colnames(mat))) {
    return(NULL)
  }
  beta <- mat[, bcols, drop = FALSE]
  eta <- beta %*% t(design$X)                          # S x N linear predictor
  family$linkinv(eta)  # S x N response-scale mean
}

#' Linear-pool envelope diagnostic (Stage-0 SS9.2)
#'
#' An equal-weight mixture of the contributing lenses' conditional-mean
#' predictives, deliberately **over-dispersed** to serve as an honest reference
#' bracket rather than a confident prediction. For each observation the envelope
#' half-width combines the average within-lens posterior uncertainty with the
#' between-lens disagreement, the latter inflated by `inflate` (default 2). The
#' reported envelope is the full equal-weight pool; flagging is leave-one-out
#' --- a lens is flagged at an observation when its posterior-mean response
#' falls outside the over-dispersed bracket formed by the *other* lenses, so a
#' disagreement carrier is caught against the consensus rather than inflating
#' its own bracket. The envelope is a diagnostic surface only and is never
#' shipped as the opinion (I2).
#'
#' Operates on the conditional-mean scale `mu_i = linkinv(eta_i)`, which is
#' comparable across every lens regardless of likelihood-tail differences, so a
#' robust-t lens and a Gaussian lens are bracketed on the same footing.
#'
#' @param fb_fits A list of contributed `koine_fb_fit` objects (opinion lenses).
#' @param pass The originating `koine_pass`.
#' @param inflate Over-dispersion multiplier on the between-lens disagreement
#'   variance (default 2, per SS9.2's "inflated approximately 2x disagreement").
#' @param level Central interval width for the bracket (default 0.95).
#' @return A list with `obs_id`, `pool_mean`, `lower`, `upper` (length-N
#'   vectors), `flags` (a data frame of `lens_id` x out-of-envelope counts),
#'   `n_components`, and `inflate`; or `NULL` when fewer than two lenses
#'   expose a comparable mean-response (no envelope is meaningful).
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = 1 + rnorm(50), x = rnorm(50))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fits <- list(
#'   fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L),
#'   fit_robust_t(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L))
#' env <- linear_pool_envelope(fits, pass)
#' env$flags
#' }
#' @export
linear_pool_envelope <- function(fb_fits, pass, inflate = 2, level = 0.95) {
  design <- .pass_design(pass)
  obs_id <- pass$obs_dim_names$obs_id
  per_lens <- list()
  for (fit in fb_fits) {
    mu <- .fit_mean_response(fit, design, pass$family)
    if (!is.null(mu)) {
      per_lens[[fit$provenance$method]] <- list(
        mean = colMeans(mu),
        var = apply(mu, 2L, stats::var)
      )
    }
  }
  K <- length(per_lens)
  if (K < 2L) {
    return(NULL)
  }
  ids <- names(per_lens)
  N <- design$N
  # K x N matrices of per-lens posterior means and within-lens variances.
  M <- t(vapply(per_lens, function(z) z$mean, numeric(N)))
  W <- t(vapply(per_lens, function(z) z$var, numeric(N)))
  pool_mean <- colMeans(M)                             # equal-weight pool mean
  within <- colMeans(W)  # mean within-lens variance
  between <- apply(M, 2L, stats::var)  # between-lens disagreement
  env_sd <- sqrt(within + inflate * between)           # over-dispersed bracket
  z <- stats::qnorm(1 - (1 - level) / 2)
  lower <- pool_mean - z * env_sd
  upper <- pool_mean + z * env_sd
  # Flag, per lens, the observations whose posterior-mean response sits outside
  # the over-dispersed bracket built from the OTHER lenses (leave-one-out). A
  # lens cannot inflate its own bracket this way, so a disagreement carrier is
  # detected against the consensus of the rest rather than hidden by it.
  outside <- vapply(seq_len(K), function(k) {
    others <- setdiff(seq_len(K), k)
    m_o <- M[others, , drop = FALSE]
    pm <- colMeans(m_o)
    wi <- colMeans(W[others, , drop = FALSE])
    bw <- if (length(others) > 1L) apply(m_o, 2L, stats::var) else 0
    sd_o <- sqrt(wi + inflate * bw)
    sum(M[k, ] < pm - z * sd_o | M[k, ] > pm + z * sd_o)
  }, integer(1))
  flags <- data.frame(lens_id = ids, n_outside = outside,
                      n_obs = N, stringsAsFactors = FALSE)
  list(obs_id = obs_id, pool_mean = pool_mean, lower = lower, upper = upper,
       flags = flags, n_components = K, inflate = inflate, level = level)
}

#' Evidence-axis diagnostic (Stage-0 SS9.3)
#'
#' Per-lens log marginal likelihood (model evidence) via bridge sampling
#' (Meng & Wong 1996; Gronau et al. 2017), as a Bayes-factor cross-check on the
#' battery. This axis has `role = "diagnostic"`: it is a scalar per lens and
#' **never enters stacking** (the synthesis assembly excludes diagnostic-role
#' output). The estimate is **prior-conditional** --- it is computed under
#' koine's canonical reference prior (coefficients N(0, `beta_prior_sd`);
#' Gaussian dispersion via a log-sigma N(0, `log_sigma_prior_sd`) prior), and
#' Bayes factors are sensitive to that prior, so the result is flagged
#' accordingly. Lenses that fit an alternative model structure (a different
#' likelihood, e.g. robust-t) are excluded from the cross-check: their evidence
#' is not comparable on the shared prior.
#'
#' Requires the suggested \pkg{bridgesampling} package; absent it, returns a
#' structured unavailable record rather than failing (I1).
#'
#' @param fb_fits A list of contributed `koine_fb_fit` objects.
#' @param pass The originating `koine_pass`.
#' @param beta_prior_sd Prior SD on the coefficients (default 5, matching the
#'   lenses' weakly-informative default).
#' @param log_sigma_prior_sd Prior SD on log-sigma for Gaussian families
#'   (default 2).
#' @param log_tau_prior_sd Prior SD on log-tau for mixed-model passes (default
#'   2, matching the RE lenses' weakly-informative default).
#' @return A list with `available` (logical), `prior_conditional` (always TRUE),
#'   and, when available, `evidence` (a data frame of `lens_id`, `log_ml`,
#'   `log_ml_se`) plus `note`. `role` is documented as `"diagnostic"`.
#' @examplesIf requireNamespace("bridgesampling", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rnorm(40), x = rnorm(40))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit <- fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
#' evidence_axis(list(fit), pass)$evidence
#' }
#' @export
evidence_axis <- function(fb_fits, pass, beta_prior_sd = 5,
                          log_sigma_prior_sd = 2, log_tau_prior_sd = 2) {
  if (!requireNamespace("bridgesampling", quietly = TRUE)) {
    return(list(available = FALSE, prior_conditional = TRUE,
                role = "diagnostic",
                note = "bridgesampling not installed; evidence axis skipped."))
  }
  if (isTRUE(pass$has_re) && identical(pass$re_structure, "met")) {
    # The single-variance-component bridge marginal below cannot score a MET
    # model's two crossed/nested variance components; the diagnostic declines the
    # MET case honestly rather than score the wrong (single-RE) marginal (I1).
    return(list(available = FALSE, prior_conditional = TRUE,
                role = "diagnostic",
                note = paste0("evidence axis not available for the MET ",
                              "(two crossed/nested intercepts) structure.")))
  }
  design <- .pass_design(pass)
  has_sigma <- identical(pass$family$name, "gaussian")
  bcols <- paste0("b_", colnames(design$X))
  # For a mixed-model pass the comparable evidence is the RE MARGINAL
  # likelihood (the RE integrated out under each lens's tau), the quantity both
  # RE lenses' (beta, sigma, tau) draws actually target --- never the population
  # likelihood, which would score the wrong model (I1).
  re <- isTRUE(pass$has_re)
  redesign <- if (re) .pass_re_design(pass) else NULL
  gh <- if (re) .gh_rule() else NULL

  log_post_unnorm <- function(theta, data) {
    beta <- theta[bcols]
    eta <- as.numeric(design$X %*% beta)
    aux <- list()
    lp_extra <- 0
    if (has_sigma) {
      log_sigma <- theta[["log_sigma"]]
      aux$sigma <- exp(log_sigma)
      lp_extra <- stats::dnorm(log_sigma, 0, log_sigma_prior_sd, log = TRUE)
    }
    if (re) {
      log_tau <- theta[["log_tau"]]
      lp_extra <- lp_extra +
        stats::dnorm(log_tau, 0, log_tau_prior_sd, log = TRUE)
      ll <- .re_group_marginal_loglik(pass$family, design$y, eta,
                                      redesign$group, exp(log_tau), aux, gh)
    } else {
      ll <- sum(pass$family$loglik_pointwise(design$y, eta, aux))
    }
    prior_beta <- sum(stats::dnorm(beta, 0, beta_prior_sd, log = TRUE))
    ll + prior_beta + lp_extra
  }

  rows <- list()
  for (fit in fb_fits) {
    if (identical(fit$status, "abstained") ||
          identical(fit$provenance$role, "diagnostic")) {
      next
    }
    # Only lenses fitting the pass family on the shared prior are comparable;
    # structure-variant lenses (own likelihood) are excluded from the axis (I1).
    if (isTRUE(fit$provenance$structure_variant)) {
      next
    }
    mat <- unclass(posterior::as_draws_matrix(fit$draws))
    if (!all(bcols %in% colnames(mat))) {
      next
    }
    lb <- stats::setNames(rep(-Inf, length(bcols)), bcols)
    ub <- stats::setNames(rep(Inf, length(bcols)), bcols)
    samples <- mat[, bcols, drop = FALSE]
    if (has_sigma) {
      if (!"sigma" %in% colnames(mat)) next
      samples <- cbind(samples, log_sigma = log(mat[, "sigma"]))
      lb <- c(lb, log_sigma = -Inf)
      ub <- c(ub, log_sigma = Inf)
    }
    if (re) {
      if (!"tau" %in% colnames(mat)) next
      samples <- cbind(samples, log_tau = log(mat[, "tau"]))
      lb <- c(lb, log_tau = -Inf)
      ub <- c(ub, log_tau = Inf)
    }
    br <- tryCatch(
      bridgesampling::bridge_sampler(
        samples = samples, log_posterior = log_post_unnorm,
        data = list(), lb = lb, ub = ub, silent = TRUE),
      error = function(e) NULL)
    if (is.null(br)) next
    rows[[length(rows) + 1L]] <- data.frame(
      lens_id = fit$provenance$method,
      log_ml = as.numeric(bridgesampling::logml(br)),
      log_ml_se = tryCatch(as.numeric(bridgesampling::error_measures(br)$re2),
                           error = function(e) NA_real_),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(list(available = FALSE, prior_conditional = TRUE,
                role = "diagnostic",
                note = "no lens yielded a comparable evidence estimate."))
  }
  list(available = TRUE, prior_conditional = TRUE, role = "diagnostic",
       evidence = do.call(rbind, rows),
       note = paste("Bayes factors are prior-conditional (beta ~ N(0,",
                    beta_prior_sd, ")); structure-variant lenses excluded."))
}
