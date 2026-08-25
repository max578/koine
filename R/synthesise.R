## ----------------------------------------------------------------------------
## The synthesis contract (Stage-0 SS6) + the two-tier gate (SS7)
##
## Predictive-space LOO-stacking over opinion lenses, PSIS-gated. Combination
## happens on the predictive density (I2); parameter posteriors are never
## averaged across heterogeneous structures. Diagnostic-role lenses are excluded
## from assembly (they never enter stacking).
## ----------------------------------------------------------------------------

#' Assemble per-lens pointwise log-likelihood matrices (Stage-0 SS6.1)
#'
#' For each contributed opinion fit, takes the native pointwise log-likelihood
#' if present, else reconstructs the S x N matrix by evaluating the family's
#' `loglik_pointwise` at each draw's linear predictor. Columns are reindexed to
#' the pass's canonical `obs_id` order (A3).
#'
#' @param fb_fits A list of `koine_fb_fit` objects.
#' @param pass The originating `koine_pass`.
#' @return A named list (by lens method) of `list(ll =, chain_id =)`, where `ll`
#'   is an S x N matrix over the canonical observation axis.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rnorm(40), x = rnorm(40))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
#' ll <- assemble_loglik(list(fit), pass)
#' dim(ll[[1]]$ll)
#' }
#' @export
assemble_loglik <- function(fb_fits, pass) {
  design <- .pass_design(pass)
  obs_id <- pass$obs_dim_names$obs_id
  out <- list()
  for (fit in fb_fits) {
    if (identical(fit$status, "abstained") ||
          identical(fit$provenance$role, "diagnostic")) {
      next
    }
    # If a lens emits its OWN pointwise log-likelihood (a structure-variant
    # lens computes its predictive under its own model, e.g. a Student-t
    # likelihood), use it directly (SS6.1) rather than the pass family.
    if (!is.null(fit$log_likelihood)) {
      out[[fit$provenance$method]] <- .native_loglik_to_list(
        fit$log_likelihood, obs_id)
      next
    }
    df <- posterior::as_draws_df(fit$draws)
    chain_id <- as.integer(df$.chain)
    mat <- unclass(posterior::as_draws_matrix(fit$draws))
    bcols <- paste0("b_", colnames(design$X))
    if (!all(bcols %in% colnames(mat))) {
      koine_abort("obs_axis_mismatch",
                  "draws lack the expected b_<col> coefficient variables.")
    }
    beta <- mat[, bcols, drop = FALSE]
    eta <- beta %*% t(design$X)                       # S x N fixed predictor
    has_sigma <- "sigma" %in% colnames(mat)
    # Random-effects branch: stack on the per-observation INTEGRATED predictive
    # (the random effect marginalised under each draw's variance components), not
    # the conditional log-likelihood --- so every RE lens predicts the identical
    # held-out quantity and the PSIS gate stays in its reliable regime (Stage-2).
    re <- isTRUE(pass$has_re)
    # MET draws carry one `tau_<grouping>` column per crossed/nested intercept;
    # the single-RE case carries a lone `tau`. The integrated predictive for a new
    # observation marginalises a fresh effect from every random factor (Stage-3).
    tau_cols <- grep("^tau", colnames(mat), value = TRUE)
    met <- re && length(tau_cols) > 1L
    if (re && !met && !("tau" %in% colnames(mat))) {
      koine_abort("obs_axis_mismatch",
                  "random-effects draws lack the expected tau variable.")
    }
    gh <- if (re && !met) .gh_rule() else NULL
    tau <- if (re && !met) mat[, "tau"] else NULL
    tau_mat <- if (met) mat[, tau_cols, drop = FALSE] else NULL
    S <- nrow(mat)
    ll <- matrix(NA_real_, nrow = S, ncol = design$N)
    for (s in seq_len(S)) {
      aux <- if (has_sigma) list(sigma = mat[s, "sigma"]) else list()
      ll[s, ] <- if (met) {
        .re_integrated_loglik_met(design$y, eta[s, ], tau_mat[s, ],
                                  mat[s, "sigma"])
      } else if (re) {
        .re_integrated_loglik(pass$family, design$y, eta[s, ], tau[s], aux, gh)
      } else {
        pass$family$loglik_pointwise(design$y, eta[s, ], aux)
      }
    }
    colnames(ll) <- obs_id                             # A3 alignment
    out[[fit$provenance$method]] <- list(ll = ll, chain_id = chain_id)
  }
  out
}

# Convert a native log_likelihood array[chain, draw, obs] into the
# list(ll = S x N matrix, chain_id) shape the gate consumes, chain-major.
.native_loglik_to_list <- function(A, obs_id) {
  dd <- dim(A)
  C <- dd[1L]
  D <- dd[2L]
  Nobs <- dd[3L]
  ll <- matrix(aperm(A, c(2L, 1L, 3L)), nrow = C * D, ncol = Nobs)
  colnames(ll) <- obs_id
  list(ll = ll, chain_id = rep(seq_len(C), each = D))
}

#' The two-tier gate (Stage-0 SS7) --- tier 1 (PSIS reliability)
#'
#' Computes per-lens PSIS-LOO with a sample-size-dependent Pareto-k threshold
#' `min(1 - 1 / log10(S), 0.7)`. Lenses whose worst k exceeds the threshold are
#' dropped with a recorded reason (never silently). This build implements the
#' drop action; K-fold refit is a later stage.
#'
#' @param loglik_list Output of [assemble_loglik()].
#' @return A list with `loo` (named list of loo objects for retained lenses),
#'   `kept` (ids), `dropped` (`list(id =, reason =, k =)`), and `threshold`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rnorm(40), x = rnorm(40))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
#' gate_psis(assemble_loglik(list(fit), pass))
#' }
#' @export
gate_psis <- function(loglik_list) {
  kept <- list()
  dropped <- list()
  threshold <- NA_real_
  for (id in names(loglik_list)) {
    ll <- loglik_list[[id]]$ll
    chain_id <- loglik_list[[id]]$chain_id
    S <- nrow(ll)
    threshold <- min(1 - 1 / log10(S), 0.7)
    r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id)
    lo <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
    max_k <- max(loo::pareto_k_values(lo))
    if (max_k >= threshold) {
      dropped[[length(dropped) + 1L]] <-
        list(id = id, reason = "psis_unreliable", k = max_k)
    } else {
      kept[[id]] <- lo
    }
  }
  list(loo = kept, kept = names(kept), dropped = dropped,
       threshold = threshold)
}

#' LOO-stacking of opinion lenses (Stage-0 SS6.2)
#'
#' Computes stacking weights over the retained lenses and the stacked
#' leave-one-out predictive density (the stacking objective evaluated at the
#' optimum). Weights are retained for the diagnostic surface; they are never
#' shipped as a pooled posterior (I2).
#'
#' @param loo_list A named list of loo objects (the gate's `loo`).
#' @return A list with `weights` (named simplex), `elpd_stacked`, and
#'   `n_components`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = 1 + rnorm(60), x = rnorm(60))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fits <- list(
#'   fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L),
#'   fit_robust_t(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L))
#' gate <- gate_psis(assemble_loglik(fits, pass))
#' if (length(gate$kept) >= 2L) stack_opinion(gate$loo)
#' }
#' @export
stack_opinion <- function(loo_list) {
  if (length(loo_list) < 2L) {
    stop("stack_opinion(): need >= 2 retained lenses to stack.", call. = FALSE)
  }
  # N x K matrix of LOO log predictive densities (elpd_loo pointwise) per lens.
  lpd <- vapply(loo_list, function(l) l$pointwise[, "elpd_loo"],
                numeric(nrow(loo_list[[1L]]$pointwise)))
  # LOO-stacking of predictive distributions (Yao et al. 2018). Optimised
  # directly rather than via loo::loo_model_weights(), because that routine
  # requires an identical posterior draw count across models --- which a
  # heterogeneous battery (different samplers, different chain lengths) violates
  # by design. The objective maximised here is identical: the weights w on the
  # simplex maximising sum_i log sum_k w_k exp(lpd_{k,i}).
  w <- stats::setNames(.stacking_weights(lpd), names(loo_list))
  mixed <- apply(lpd, 1L, function(row) .log_sum_exp(row + log(w)))
  list(weights = w, elpd_stacked = sum(mixed), n_components = length(loo_list))
}

.log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.softmax <- function(z) {
  z <- z - max(z)
  e <- exp(z)
  e / sum(e)
}

# Stacking weights: maximise the LOO log predictive density of the mixture over
# the simplex. Reparameterised by softmax (first component fixed) so the
# optimisation is unconstrained; K is small (a handful of lenses).
.stacking_weights <- function(lpd) {
  K <- ncol(lpd)
  if (K == 1L) return(1)
  lpd_s <- lpd - apply(lpd, 1L, max)        # per-obs stabilisation (shift)
  P <- exp(lpd_s)                            # N x K predictive densities
  neg_obj <- function(z) {
    w <- .softmax(c(0, z))
    -sum(log(pmax(as.numeric(P %*% w), 1e-300)))
  }
  opt <- stats::optim(rep(0, K - 1L), neg_obj, method = "BFGS")
  .softmax(c(0, opt$par))
}

#' Weight-concentration diagnostic (Stage-0 SS9)
#'
#' Summarises how concentrated the stacking weights are. Concentration on one
#' component means the battery added no predictive diversity (the lenses
#' corroborate rather than diversify); spread weights mean genuine
#' model-structure diversity was realised. Reports the effective number of
#' components `n_eff = 1 / sum(w^2)`, the maximum weight, the entropy, and a
#' `degenerate` flag.
#'
#' @param weights A named numeric vector of stacking weights.
#' @param n_eff_threshold Below this effective component count the weights are
#'   flagged degenerate. Default 1.5.
#' @return A list with `weights`, `max_weight`, `n_eff`, `entropy`, and
#'   `degenerate`.
#' @examples
#' # spread weights: genuine diversity; concentrated weights: corroboration
#' weight_concentration(c(dram = 0.55, robust_t = 0.45))
#' weight_concentration(c(dram = 0.97, robust_t = 0.03))
#' @export
weight_concentration <- function(weights, n_eff_threshold = 1.5) {
  w <- weights / sum(weights)
  n_eff <- 1 / sum(w^2)
  ent <- -sum(ifelse(w > 0, w * log(w), 0))
  list(weights = w, max_weight = max(w), n_eff = n_eff, entropy = ent,
       degenerate = n_eff < n_eff_threshold)
}

#' Synthesise the fourth opinion for a pass (Stage-0 SS1 flow)
#'
#' Resolves the battery, fits each opinion lens, assembles pointwise
#' log-likelihoods, applies the PSIS gate, stacks the survivors, and reports a
#' weight-concentration diagnostic --- or returns a first-class [abstention()]
#' when fewer than two independent opinion lenses survive. The result is either
#' a `koine_result` or a `koine_abstention`, never a silent partial.
#'
#' @param pass A `koine_pass`.
#' @param registry A registry list; defaults to [koine_stage1_registry()].
#' @param on_degenerate What to do when the stacking weights are degenerate
#'   (the battery corroborated rather than diversified): `"flag"` (default)
#'   returns the opinion with the degeneracy recorded in its diagnostics;
#'   `"abstain"` returns a `weight_degenerate` abstention (SS8).
#' @param keep_fits If `TRUE`, attach the contributing per-lens fits to the
#'   result (its `fits` field), so [koine_opinion_draws()] can materialise the
#'   stacking-mixture posterior without re-fitting. Default `FALSE` (the
#'   diagnostic-only contract is unchanged).
#' @return A `koine_result` or a `koine_abstention`.
#' @export
#' @examples
#' d <- data.frame(y = rnorm(30), x = rnorm(30))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' # with a single-lens battery, koine abstains rather than guessing
#' res <- synthesise(pass, registry = list(koine_stage1_registry()[[1]]))
#' is_abstention(res)
synthesise <- function(pass, registry = koine_stage1_registry(),
                       on_degenerate = c("flag", "abstain"),
                       keep_fits = FALSE) {
  on_degenerate <- match.arg(on_degenerate)
  stopifnot(inherits(pass, "koine_pass"))
  if (isTRUE(pass$has_re) && !isTRUE(pass$re_supported)) {
    return(abstention(
      "feature_unsupported", scope = "synthesis",
      detail = paste0(
        "the random-effects structure is outside this build's support ",
        "(a single random intercept on one grouping factor); koine abstains ",
        "rather than mis-modelling it.")))
  }
  resolved <- resolve_battery(pass, registry)
  opinion <- Filter(function(L) identical(L$role, "opinion"),
                    resolved$enrolled)
  if (length(opinion) < 2L) {
    return(abstention(
      "no_diverse_battery", scope = "synthesis",
      detail = sprintf(
        "only %d opinion lens(es) resolved; >= 2 independent lenses needed.",
        length(opinion)),
      lenses_dropped = resolved$excluded))
  }
  fits <- Map(function(L, fit) list(lens = L, fit = fit), opinion,
             lapply(opinion, function(L) L$fit_fn(pass, L)))
  # A lens declining at fit time returns a `koine_abstention` value (never a
  # condition, never a bare `NA`) --- carry its true reason forward rather than
  # silently dropping it and letting the PSIS gate below take the blame (K1).
  lens_abstained <- Filter(Negate(is.null), lapply(fits, function(pair) {
    if (is_abstention(pair$fit)) {
      list(id = pair$lens$id, reason = pair$fit$reason)
    } else {
      NULL
    }
  }))
  contributed <- Filter(
    function(pair) identical(pair$fit$status, "contributed"), fits)
  fits <- lapply(contributed, `[[`, "fit")
  if (length(fits) < 2L) {
    reasons <- unique(vapply(lens_abstained, `[[`, character(1L), "reason"))
    # A single distinct lens reason propagates as the synthesis-level reason;
    # a mix of reasons (or none, e.g. a lens contributed but was later dropped
    # for some other cause) falls back to the generic diversity abstention.
    reason <- if (length(reasons) == 1L) reasons else "no_diverse_battery"
    return(abstention(
      reason, scope = "synthesis",
      detail = sprintf(
        "only %d of %d opinion lens(es) contributed a fit; >= 2 needed.",
        length(fits), length(opinion)),
      lenses_dropped = c(resolved$excluded, lens_abstained)))
  }
  loglik_list <- assemble_loglik(fits, pass)
  gate <- gate_psis(loglik_list)
  if (length(gate$kept) < 2L) {
    return(abstention(
      "psis_unreliable", scope = "synthesis",
      detail = "fewer than two lenses survived the PSIS reliability gate.",
      lenses_dropped = c(lens_abstained, gate$dropped)))
  }
  stacked <- stack_opinion(gate$loo)
  conc <- weight_concentration(stacked$weights)
  if (conc$degenerate && identical(on_degenerate, "abstain")) {
    return(abstention(
      "weight_degenerate", scope = "synthesis",
      detail = sprintf(
        "only %.2f of %d effective lens(es) carry weight; diversity lost.",
        conc$n_eff, stacked$n_components)))
  }
  # The diagnostic layer (SS9): the in-stack concentration check, the
  # over-dispersed linear-pool envelope, and the prior-conditional evidence
  # axis. All three are surfaced, none is ever shipped as the opinion (I2).
  diagnostics <- list(
    weight_concentration = conc,
    linear_pool_envelope = linear_pool_envelope(fits, pass),
    evidence_axis = evidence_axis(fits, pass)
  )
  structure(
    list(opinion = stacked, gate = gate, diagnostics = diagnostics,
         independence = independence_report(fits),
         excluded = resolved$excluded, ir_hash = pass$ir_hash,
         fits = if (isTRUE(keep_fits)) fits else NULL),
    class = "koine_result"
  )
}

#' @export
print.koine_result <- function(x, ...) {
  cat("<koine result --- synthesised fourth opinion>\n")
  cat("  components: ", x$opinion$n_components, "\n", sep = "")
  cat("  weights   :\n")
  for (nm in names(x$opinion$weights)) {
    cat(sprintf("    %-14s %.3f\n", nm, x$opinion$weights[[nm]]))
  }
  cat("  elpd (stacked LOO): ", round(x$opinion$elpd_stacked, 2), "\n",
      sep = "")
  cat("  PSIS k threshold  : ", round(x$gate$threshold, 3), "\n", sep = "")
  wc <- x$diagnostics$weight_concentration
  cat("  diversity         : ", round(wc$n_eff, 2), " effective lens(es)",
      if (wc$degenerate) "  [DEGENERATE: corroboration, not diversity]" else "",
      "\n", sep = "")
  lpe <- x$diagnostics$linear_pool_envelope
  if (!is.null(lpe)) {
    flagged <- sum(lpe$flags$n_outside > 0L)
    cat("  pool envelope     : ", lpe$n_components, " lenses, ", flagged,
        " disagreeing on >=1 obs (", lpe$level * 100, "%, ", lpe$inflate,
        "x inflated)\n", sep = "")
  }
  ev <- x$diagnostics$evidence_axis
  if (isTRUE(ev$available)) {
    cat("  evidence axis     : ", nrow(ev$evidence),
        " lens(es), log-ml [prior-conditional]\n", sep = "")
  }
  invisible(x)
}
