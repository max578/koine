## ----------------------------------------------------------------------------
## The validation / maturity-gate contract (Stage-0 SS11, open-arms SS6)
##
## A koine-first lens vouches for its own correctness, which weakens the
## "agreement across independent implementations is evidence" story. So a
## koine-authored lens must clear an assessment gate before it may contribute to
## triangulation, and it carries an honest maturity tier:
##
##   reference-validated  >  koine-implemented-validated  >  experimental
##
## Promotion of an experimental lens to koine-implemented-validated requires
## (SS11): simulation-based calibration (SBC), a reference-posterior agreement
## check, and cross-engine concordance. The first two are runnable in this
## standalone, pre-graduation package; cross-engine concordance needs the engine
## panel and is an explicit scope boundary deferred to graduation (recorded,
## never silently skipped). A lens lacking validation evidence is treated as
## experimental by the gate and abstains --- it is never silently trusted.
## ----------------------------------------------------------------------------

#' Simulation-based calibration of a lens (Stage-0 SS11; Talts et al. 2018)
#'
#' Repeatedly draws true parameters from the prior, simulates a dataset from the
#' lens's own generative model, fits the lens, and records the rank of each true
#' parameter among the (thinned) posterior draws. If the lens samples its
#' posterior correctly, the ranks are uniform on `0:n_keep`; systematic
#' departures (a chi-squared test per parameter) reveal a miscalibrated sampler.
#'
#' SBC validates the `{prior, model, sampler}` triple as a whole, so the
#' `generate` callback must simulate data from exactly the model the lens
#' assumes (e.g. a Student-t likelihood for the robust-t lens), using the
#' exact prior the lens places.
#'
#' @param generate A function `(seed) -> list(theta = named-numeric, pass =
#'   koine_pass)`: draws true parameters from the prior and builds a pass whose
#'   data is simulated from the model at `theta`. `theta`'s names must match the
#'   posterior draw variables to be ranked.
#' @param fit_fn A function `(pass) -> koine_fb_fit` (the lens, with its sampler
#'   controls already set).
#' @param param_names Character vector of parameters to rank (a subset of the
#'   names of `theta` and of the fit's draw variables).
#' @param n_sims Number of SBC replications.
#' @param n_keep Posterior draws retained per fit (the ranks lie in `0:n_keep`);
#'   draws are thinned to this many to keep the rank statistic near-independent.
#' @param n_bins Number of histogram bins for the chi-squared uniformity test.
#' @param alpha Family-wise significance level (Bonferroni-adjusted across
#'   parameters) below which calibration is judged to have failed.
#' @param seed Master seed; per-replication seeds are derived from it.
#' @return A list with `ranks` (an `n_sims` x `length(param_names)` matrix),
#'   `n_keep`, `uniformity` (a data frame of `param`, `chisq`, `df`, `p_value`),
#'   `alpha_adj`, and `passed` (logical: every parameter's p-value exceeds the
#'   adjusted alpha).
#' @examples
#' \donttest{
#' setup <- koine_lens_validation_setup(
#'   "dram", n_obs = 30L,
#'   fit_control = list(n_iter = 800L, n_burn = 400L, n_chains = 2L))
#' res <- sbc_check(setup$generate, setup$fit_fn, setup$param_names,
#'                  n_sims = 16L, n_keep = 16L, n_bins = 4L)
#' res$passed
#' }
#' @export
sbc_check <- function(generate, fit_fn, param_names, n_sims = 128L,
                      n_keep = 64L, n_bins = 8L, alpha = 0.05, seed = 1L) {
  ranks <- matrix(NA_real_, nrow = n_sims, ncol = length(param_names),
                  dimnames = list(NULL, param_names))
  for (s in seq_len(n_sims)) {
    g <- generate(seed + s)
    fit <- fit_fn(g$pass)
    if (is_abstention(fit) || !identical(fit$status, "contributed")) {
      next
    }
    mat <- unclass(posterior::as_draws_matrix(fit$draws))
    idx <- unique(round(seq(1, nrow(mat), length.out = n_keep)))
    thin <- mat[idx, param_names, drop = FALSE]
    for (p in param_names) {
      ranks[s, p] <- sum(thin[, p] < g$theta[[p]])
    }
  }
  ranks <- ranks[stats::complete.cases(ranks), , drop = FALSE]
  unif <- .sbc_uniformity(ranks, n_keep, n_bins)
  alpha_adj <- alpha / length(param_names)
  list(ranks = ranks, n_keep = n_keep, uniformity = unif, alpha_adj = alpha_adj,
       n_eff_sims = nrow(ranks),
       passed = all(unif$p_value > alpha_adj))
}

# Per-parameter chi-squared test that SBC ranks are uniform on 0:n_keep.
.sbc_uniformity <- function(ranks, n_keep, n_bins) {
  rows <- lapply(colnames(ranks), function(p) {
    r <- ranks[, p]
    br <- seq(0, n_keep + 1, length.out = n_bins + 1L)
    obs <- as.numeric(table(cut(r, breaks = br, include.lowest = TRUE,
                                right = FALSE)))
    expct <- rep(length(r) / n_bins, n_bins)
    chisq <- sum((obs - expct)^2 / expct)
    df <- n_bins - 1L
    data.frame(param = p, chisq = chisq, df = df,
               p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Exact grid-based reference-posterior check (Stage-0 SS11)
#'
#' Compares a lens's posterior to an exact reference computed by dense grid
#' integration of the same unnormalised log-posterior --- a reference that is
#' independent of any sampler and needs no external engine, so it is runnable in
#' this standalone package. Feasible only for low-dimensional models (the grid
#' is the product over parameters); for larger models it returns a structured
#' "not applicable" record rather than a misleading number (I1).
#'
#' @param fit A contributed `koine_fb_fit` from the lens under test.
#' @param log_post_unnorm A function `(theta) -> numeric`: the lens's
#'   unnormalised log-posterior (likelihood + prior) at a named parameter vector
#'   on the **draw scale** (the same variable names the fit's draws carry).
#' @param grid A named list of numeric grid vectors, one per parameter, naming
#'   the parameters in draw-scale order.
#' @param tol Maximum tolerated absolute standardised discrepancy between the
#'   lens and reference posterior mean (in reference SDs) for a pass.
#' @param max_cells Refuse (return not-applicable) above this many grid cells.
#' @return A list with `applicable` (logical), and when applicable `reference`
#'   and `lens` summaries (per-parameter mean and SD), `discrepancy` (per
#'   parameter, in reference SDs), and `passed`.
#' @examples
#' \donttest{
#' setup <- koine_lens_validation_setup(
#'   "dram", n_obs = 30L,
#'   fit_control = list(n_iter = 1500L, n_burn = 750L, n_chains = 2L))
#' g <- setup$generate(1L)
#' fit <- setup$fit_fn(g$pass)
#' ref <- setup$reference(g, n_pts = 21L)
#' reference_posterior_check(fit, ref$log_post_unnorm, ref$grid)$passed
#' }
#' @export
reference_posterior_check <- function(fit, log_post_unnorm, grid, tol = 0.15,
                                      max_cells = 2e5) {
  pnames <- names(grid)
  n_cells <- prod(vapply(grid, length, integer(1)))
  if (n_cells > max_cells) {
    return(list(applicable = FALSE,
                note = sprintf(
                  "model too large for a grid reference (%g cells > %g).",
                  n_cells, max_cells)))
  }
  cells <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE)
  lp <- apply(cells, 1L, function(row) {
    log_post_unnorm(stats::setNames(as.numeric(row), pnames))
  })
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  ref_mean <- vapply(pnames, function(p) sum(w * cells[[p]]), numeric(1))
  ref_sd <- vapply(pnames, function(p) {
    sqrt(sum(w * (cells[[p]] - ref_mean[[p]])^2))
  }, numeric(1))
  mat <- unclass(posterior::as_draws_matrix(fit$draws))
  lens_mean <- vapply(pnames, function(p) mean(mat[, p]), numeric(1))
  lens_sd <- vapply(pnames, function(p) stats::sd(mat[, p]), numeric(1))
  discrepancy <- abs(lens_mean - ref_mean) / ref_sd
  list(applicable = TRUE,
       reference = data.frame(param = pnames, mean = ref_mean, sd = ref_sd,
                              stringsAsFactors = FALSE),
       lens = data.frame(param = pnames, mean = lens_mean, sd = lens_sd,
                         stringsAsFactors = FALSE),
       discrepancy = stats::setNames(discrepancy, pnames),
       tol = tol, passed = all(discrepancy < tol))
}

#' The corroboration / maturity gate (Stage-0 SS11)
#'
#' Decides whether a lens of a given maturity tier may corroborate the fourth
#' opinion, and at what weight discount. A koine-first lens
#' (`koine-implemented-validated`) may corroborate only when validation evidence
#' is on record (`validation_ref`); without it the lens is treated as
#' `experimental` and may not corroborate --- an unvalidated first-mover is
#' never silently trusted (open-arms correctness, SS6/SS11).
#'
#' @param maturity_tier One of `"reference-validated"`,
#'   `"koine-implemented-validated"`, `"experimental"`.
#' @param validation_ref Optional identifier of the validation evidence (e.g. a
#'   ledger cairn id). Required to back the `koine-implemented-validated` tier.
#' @return A list with `effective_tier`, `may_corroborate` (logical), and
#'   `discount` (1 reference-validated, 0.5 koine-validated, 0 experimental).
#' @examples
#' gate_corroboration("reference-validated")
#' gate_corroboration("koine-implemented-validated", validation_ref = "cairn-1")
#' gate_corroboration("koine-implemented-validated") # unbacked -> experimental
#' @export
gate_corroboration <- function(maturity_tier, validation_ref = NULL) {
  effective <- maturity_tier
  if (identical(maturity_tier, "koine-implemented-validated") &&
        (is.null(validation_ref) || !nzchar(as.character(validation_ref)))) {
    effective <- "experimental"
  }
  discount <- switch(effective,
                     "reference-validated" = 1,
                     "koine-implemented-validated" = 0.5,
                     "experimental" = 0,
                     0)
  list(effective_tier = effective,
       may_corroborate = discount > 0,
       discount = discount)
}

#' Effective maturity tier of a lens after the validation gate
#'
#' Resolves a lens spec's declared maturity against its `validation_ref`: a
#' koine-first lens without recorded validation evidence is demoted to
#' `experimental` (it may not corroborate). Reference-validated lenses and
#' validation-backed koine lenses keep their declared tier.
#'
#' @param lens_spec A `koine_lens_spec`.
#' @return The effective maturity tier (a string).
#' @examples
#' lens_effective_maturity(koine_stage1_registry()[[1]])
#' @export
lens_effective_maturity <- function(lens_spec) {
  gate_corroboration(lens_spec$maturity_tier,
                     lens_spec$validation_ref)$effective_tier
}

# Whether a lens's maturity must be backed by validation evidence to corroborate
# (koine-authored tiers), as opposed to reference implementations.
.lens_needs_validation <- function(lens_spec) {
  lens_spec$maturity_tier %in%
    c("koine-implemented-validated", "experimental")
}

## -- the canonical validation harness for koine's own lenses -----------------
##
## Each koine-authored lens is validated against its OWN generative model and a
## grid-exact reference posterior. The setup below is the single source of truth
## for both the test suite (tiny replications) and the recorded validation study
## (inst/validation/sbc_study.R). It encodes, per lens: the prior, the model
## simulator, the sampler call, and the unnormalised log-posterior on the draw
## scale (including the log-sigma -> sigma Jacobian where the lens reports sigma
## but priors on log-sigma).

# Coefficient names for a K-column GLM design (intercept + x predictors).
.val_bcols <- function(n_pred) {
  c("b_(Intercept)", paste0("b_x", seq_len(n_pred)))
}

.val_formula <- function(n_pred) {
  stats::reformulate(paste0("x", seq_len(n_pred)), response = "y")
}

.val_design_X <- function(n_obs, n_pred) {
  cols <- lapply(seq_len(n_pred), function(j) stats::rnorm(n_obs))
  names(cols) <- paste0("x", seq_len(n_pred))
  list(df_pred = as.data.frame(cols),
       X = cbind(`(Intercept)` = 1, do.call(cbind, cols)))
}

# Grid centred on a quick MLE of one simulated dataset, +/- `span` SEs, with
# `n_pts` points per dimension; sigma kept positive. The window only has to
# bracket the posterior mass --- the reference integrates the exact posterior
# within it --- so coefficient SEs use the unit-information matrix with a robust
# fallback, and the expected small-sample glm.fit warnings are suppressed.
.val_grid <- function(pass, has_sigma, n_pred, n_pts, span = 4) {
  design <- .pass_design(pass)
  bcols <- .val_bcols(n_pred)
  fit <- suppressWarnings(tryCatch(
    stats::glm.fit(design$X, design$y, family = .pass_glm_family(pass)),
    error = function(e) NULL))
  if (is.null(fit)) {
    cf <- rep.int(0, design$K)
    se <- rep.int(1, design$K)
  } else {
    cf <- stats::coef(fit)
    cf[!is.finite(cf)] <- 0
    # Coefficient SEs from the IRLS info matrix (X' W X)^{-1} * dispersion.
    # --- exact for the GLM, no fragile summary.glm round-trip.
    xtwx_inv <- tryCatch(
      chol2inv(chol(crossprod(design$X * sqrt(fit$weights)))),
      error = function(e) NULL)
    disp <- if (has_sigma) {
      sum((design$y - fit$fitted.values)^2) / max(1L, design$N - design$K)
    } else {
      1
    }
    se <- if (is.null(xtwx_inv)) {
      rep.int(1, design$K)
    } else {
      sqrt(pmax(diag(xtwx_inv) * disp, 1e-6))
    }
  }
  se[!is.finite(se) | se <= 0] <- 1
  grid <- lapply(seq_len(design$K), function(k) {
    seq(cf[k] - span * se[k], cf[k] + span * se[k], length.out = n_pts)
  })
  names(grid) <- bcols
  if (has_sigma) {
    rs <- stats::sd(design$y)
    grid[["sigma"]] <- seq(max(1e-3, rs * 0.3), rs * 2.5, length.out = n_pts)
  }
  grid
}

#' Build the validation setup for a koine-authored lens
#'
#' Returns the prior draw, the model simulator, the sampler call, the parameters
#' to rank, and a grid-exact reference-posterior builder for one of koine's
#' koine-authored lenses (`"dram"`, `"pg_gibbs"`, `"robust_t"`). This is the
#' shared harness behind [sbc_check()] and [reference_posterior_check()] for the
#' recorded validation study; SBC simulates from the model returned here, so the
#' calibration test sees the lens's own `{prior, model, sampler}` triple.
#'
#' @param lens_id One of `"dram"`, `"pg_gibbs"`, `"robust_t"` (population)
#'   or `"glmm_marginal"`, `"nimble_glmm"` (mixed-model lenses).
#' @param n_obs Observations per simulated dataset (population lenses).
#' @param n_pred Number of predictors beyond the intercept.
#' @param nu Student-t degrees of freedom for `"robust_t"` (the generative model
#'   and the sampler share this value).
#' @param fit_control A list of sampler controls passed to the lens fit.
#' @param family For the mixed-model lenses, the response family
#'   (`"gaussian"`, `"binomial"`, or `"poisson"`); ignored by the population
#'   lenses, which carry their own family.
#' @param n_groups Number of random-effect groups (mixed-model lenses).
#' @param n_per_group Observations per group (mixed-model lenses).
#' @return A list with `generate`, `fit_fn`, `param_names`, and `reference`.
#' @examples
#' setup <- koine_lens_validation_setup("dram")
#' setup$param_names
#' @export
koine_lens_validation_setup <- function(lens_id, n_obs = 40L, n_pred = 1L,
                                        nu = 4,
                                        fit_control = list(n_iter = 1500L,
                                                           n_burn = 750L,
                                                           n_chains = 2L),
                                        family = c("gaussian", "binomial",
                                                   "poisson"),
                                        n_groups = 12L, n_per_group = 6L) {
  bcols <- .val_bcols(n_pred)
  beta_sd <- 5
  log_sigma_sd <- 2

  if (lens_id %in% c("glmm_marginal", "nimble_glmm", "glmm_robust_re")) {
    fam <- if (lens_id == "glmm_robust_re") "gaussian" else match.arg(family)
    return(.re_validation_setup(lens_id, fam, n_pred, bcols, beta_sd,
                                log_sigma_sd, n_groups, n_per_group,
                                fit_control, nu))
  }

  if (lens_id == "pg_gibbs") {
    has_sigma <- FALSE
    pnames <- bcols
    generate <- function(seed) {
      set.seed(seed)
      dx <- .val_design_X(n_obs, n_pred)
      beta <- stats::rnorm(n_pred + 1L, 0, beta_sd)
      eta <- as.numeric(dx$X %*% beta)
      y <- stats::rbinom(n_obs, 1L, stats::plogis(eta))
      pass <- koine_pass(.val_formula(n_pred),
                         cbind(data.frame(y = y), dx$df_pred),
                         stats::binomial(), seed = as.integer(seed))
      list(theta = stats::setNames(beta, bcols), pass = pass)
    }
    fit_fn <- function(pass) {
      do.call(fit_pg_gibbs, c(list(pass = pass), fit_control))
    }
    loglik <- function(design, theta) {
      eta <- as.numeric(design$X %*% theta[bcols])
      stats::dbinom(design$y, 1L, stats::plogis(eta), log = TRUE)
    }
  } else if (lens_id == "dram") {
    has_sigma <- TRUE
    pnames <- c(bcols, "sigma")
    generate <- function(seed) {
      set.seed(seed)
      dx <- .val_design_X(n_obs, n_pred)
      beta <- stats::rnorm(n_pred + 1L, 0, beta_sd)
      sigma <- exp(stats::rnorm(1, 0, log_sigma_sd))
      eta <- as.numeric(dx$X %*% beta)
      y <- stats::rnorm(n_obs, eta, sigma)
      pass <- koine_pass(.val_formula(n_pred),
                         cbind(data.frame(y = y), dx$df_pred),
                         stats::gaussian(), seed = as.integer(seed))
      list(theta = stats::setNames(c(beta, sigma), pnames), pass = pass)
    }
    fit_fn <- function(pass) {
      do.call(fit_dram, c(list(pass = pass), fit_control))
    }
    loglik <- function(design, theta) {
      eta <- as.numeric(design$X %*% theta[bcols])
      stats::dnorm(design$y, eta, theta[["sigma"]], log = TRUE)
    }
  } else if (lens_id == "robust_t") {
    has_sigma <- TRUE
    pnames <- c(bcols, "sigma")
    generate <- function(seed) {
      set.seed(seed)
      dx <- .val_design_X(n_obs, n_pred)
      beta <- stats::rnorm(n_pred + 1L, 0, beta_sd)
      sigma <- exp(stats::rnorm(1, 0, log_sigma_sd))
      eta <- as.numeric(dx$X %*% beta)
      y <- eta + sigma * stats::rt(n_obs, df = nu)
      pass <- koine_pass(.val_formula(n_pred),
                         cbind(data.frame(y = y), dx$df_pred),
                         stats::gaussian(), seed = as.integer(seed))
      list(theta = stats::setNames(c(beta, sigma), pnames), pass = pass)
    }
    fit_fn <- function(pass) {
      do.call(fit_robust_t, c(list(pass = pass, nu = nu), fit_control))
    }
    loglik <- function(design, theta) {
      mu <- as.numeric(design$X %*% theta[bcols])
      s <- theta[["sigma"]]
      stats::dt((design$y - mu) / s, df = nu, log = TRUE) - log(s)
    }
  } else {
    stop("koine_lens_validation_setup(): unknown lens id '", lens_id, "'.",
         call. = FALSE)
  }

  # One reference builder for every lens: the unnormalised log-posterior on the
  # draw scale, with the log-sigma -> sigma Jacobian added only for the
  # sigma-bearing lenses (so the grid reference matches the lens's reported,
  # back-transformed sigma marginal).
  reference <- function(g, n_pts = 31L) {
    design <- .pass_design(g$pass)
    lp <- function(theta) {
      val <- sum(loglik(design, theta)) +
        sum(stats::dnorm(theta[bcols], 0, beta_sd, log = TRUE))
      if (has_sigma) {
        val <- val +
          stats::dnorm(log(theta[["sigma"]]), 0, log_sigma_sd, log = TRUE) -
          log(theta[["sigma"]])
      }
      val
    }
    list(log_post_unnorm = lp,
         grid = .val_grid(g$pass, has_sigma, n_pred, n_pts))
  }
  list(generate = generate, fit_fn = fit_fn, param_names = pnames,
       reference = reference)
}

## -- mixed-model validation setup --------------------------------------------
##
## SBC for the RE lenses ranks (beta, sigma, tau) against the prior draws, with
## data simulated from the full generative model (draw tau, draw b_g ~ N(0,
## tau^2), simulate y). Both RE lenses target the *marginal* posterior p(beta,
## sigma, tau | y) --- the marginal lens by quadrature/analytic integration, the
## NIMBLE lens by sampling b and reporting only (beta, sigma, tau). Because
## marginalising b is exact, rank uniformity on (beta, sigma, tau) is correct
## SBC criterion for both. The grid-exact reference is built on the same
## random-effects MARGINAL log-posterior (`.re_group_marginal_loglik` + priors,
## with the log-sigma / log-tau Jacobians), so it never scores the population
## model by mistake (I1).

.re_validation_setup <- function(lens_id, family, n_pred, bcols, beta_sd,
                                 log_sigma_sd, n_groups, n_per_group,
                                 fit_control, nu = 4) {
  has_sigma <- family == "gaussian"
  # SBC must simulate from the lens's ACTUAL (family-appropriate) prior, so the
  # {prior, model, sampler} triple is the one being validated.
  prior <- .re_prior(family)
  beta_sd <- prior$beta_sd
  log_sigma_sd <- prior$log_sigma_sd
  log_tau_sd <- prior$log_tau_sd
  pnames <- c(bcols, if (has_sigma) "sigma", "tau")
  fam_obj <- switch(family, gaussian = stats::gaussian(),
                    binomial = stats::binomial(), poisson = stats::poisson())
  # the robust-RE lens assumes Student-t random effects; its SBC must simulate
  # the random effects from t_nu (its own generative model) and the reference
  # must use the t-RE marginal --- otherwise the {prior, model, sampler} triple
  # being validated is not the lens's own.
  is_robust <- lens_id == "glmm_robust_re"
  fit_lens <- switch(lens_id,
                     glmm_marginal = fit_glmm_marginal,
                     nimble_glmm = fit_nimble_glmm,
                     glmm_robust_re = fit_glmm_robust_re)
  re_formula <- stats::as.formula(paste(
    "y ~", paste(c(paste0("x", seq_len(n_pred)), "(1 | g)"), collapse = " + ")))

  generate <- function(seed) {
    set.seed(seed)
    n_obs <- n_groups * n_per_group
    grp <- factor(rep(seq_len(n_groups), each = n_per_group))
    dx <- .val_design_X(n_obs, n_pred)
    beta <- stats::rnorm(n_pred + 1L, 0, beta_sd)
    tau <- exp(stats::rnorm(1, 0, log_tau_sd))
    b <- if (is_robust) {
      stats::rt(n_groups, df = nu) * tau
    } else {
      stats::rnorm(n_groups, 0, tau)
    }
    eta <- as.numeric(dx$X %*% beta) + b[grp]
    sigma <- if (has_sigma) exp(stats::rnorm(1, 0, log_sigma_sd)) else NA_real_
    y <- switch(family,
                gaussian = stats::rnorm(n_obs, eta, sigma),
                binomial = stats::rbinom(n_obs, 1L, stats::plogis(eta)),
                poisson  = stats::rpois(n_obs, exp(eta)))
    df <- cbind(data.frame(y = y), dx$df_pred, g = grp)
    pass <- koine_pass(re_formula, df, fam_obj, seed = as.integer(seed))
    theta <- c(beta, if (has_sigma) sigma, tau)
    list(theta = stats::setNames(theta, pnames), pass = pass)
  }

  fit_fn <- function(pass) {
    args <- c(list(pass = pass), if (is_robust) list(nu = nu), fit_control)
    do.call(fit_lens, args)
  }

  reference <- function(g, n_pts = 15L) {
    design <- .pass_design(g$pass)
    redesign <- .pass_re_design(g$pass)
    gh <- .gh_rule()
    lp <- function(theta) {
      beta <- theta[bcols]
      eta <- as.numeric(design$X %*% beta)
      aux <- list()
      jac <- stats::dnorm(log(theta[["tau"]]), 0, log_tau_sd, log = TRUE) -
        log(theta[["tau"]])
      if (has_sigma) {
        aux$sigma <- theta[["sigma"]]
        jac <- jac +
          stats::dnorm(log(theta[["sigma"]]), 0, log_sigma_sd, log = TRUE) -
          log(theta[["sigma"]])
      }
      ll <- if (is_robust) {
        .re_group_marginal_t(design$y, eta, redesign$group, theta[["tau"]],
                             theta[["sigma"]], nu)
      } else {
        .re_group_marginal_loglik(g$pass$family, design$y, eta,
                                  redesign$group, theta[["tau"]], aux, gh)
      }
      ll + sum(stats::dnorm(beta, 0, beta_sd, log = TRUE)) + jac
    }
    list(log_post_unnorm = lp,
         grid = .re_val_grid(g$pass, has_sigma, n_pred, n_pts))
  }

  list(generate = generate, fit_fn = fit_fn, param_names = pnames,
       reference = reference)
}

# Grid for the RE reference posterior, uniformly spaced (the reference
# integrator [reference_posterior_check()] assumes a uniform measure) and
# DATA-DRIVEN so it brackets the marginal posterior tightly: a coarse,
# sd(y)-scaled window leaves the sharp mixed-model marginal in a handful of
# cells, collapsing the reference SD and exploding the standardised discrepancy.
# Coefficients centre on the population GLM fit (+/- 4 clustering-inflated SEs);
# the residual SD (sigma) and random-intercept SD (tau) centre on
# method-of-moments estimates from the group-mean decomposition. The location
# parameters (beta, sigma) are then resolved well enough to check against the
# exact grid reference; the weakly-identified tau is resolved only to grid
# accuracy, with SBC ([sbc_check()]) as its authoritative validation.
.re_val_grid <- function(pass, has_sigma, n_pred, n_pts) {
  design <- .pass_design(pass)
  redesign <- .pass_re_design(pass)
  bcols <- .val_bcols(n_pred)
  fit <- suppressWarnings(tryCatch(
    stats::glm.fit(design$X, design$y, family = .pass_glm_family(pass)),
    error = function(e) NULL))
  if (is.null(fit)) {
    cf <- rep.int(0, design$K)
    se <- rep.int(1, design$K)
  } else {
    cf <- stats::coef(fit)
    cf[!is.finite(cf)] <- 0
    xtwx_inv <- tryCatch(
      chol2inv(chol(crossprod(design$X * sqrt(fit$weights)))),
      error = function(e) NULL)
    # raw residual dispersion (it absorbs the random-effect variance) gives a
    # conservatively wide coefficient SE for the clustered marginal.
    disp <- sum((design$y - fit$fitted.values)^2) / max(1L, design$N - design$K)
    se <- if (is.null(xtwx_inv)) {
      rep.int(1, design$K)
    } else {
      sqrt(pmax(diag(xtwx_inv) * disp, 1e-4))
    }
  }
  se[!is.finite(se) | se <= 0] <- 1
  grid <- lapply(seq_len(design$K), function(k) {
    seq(cf[k] - 4 * se[k], cf[k] + 4 * se[k], length.out = n_pts)
  })
  names(grid) <- bcols
  # method-of-moments variance components from the group-mean decomposition.
  gm <- tapply(design$y - as.numeric(design$X %*% cf), redesign$group, mean)
  tau_hat <- max(0.1, stats::sd(gm))
  grid[["tau"]] <- seq(max(1e-2, tau_hat * 0.2), tau_hat * 2.6,
                       length.out = n_pts)
  if (has_sigma) {
    resid <- design$y - as.numeric(design$X %*% cf) - gm[redesign$group]
    sig_hat <- max(0.1, stats::sd(resid))
    grid[["sigma"]] <- seq(max(1e-2, sig_hat * 0.4), sig_hat * 1.9,
                           length.out = n_pts)
  }
  grid
}
