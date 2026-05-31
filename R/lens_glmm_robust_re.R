## ----------------------------------------------------------------------------
## glmm_robust_re lens --- a heavy-tailed RANDOM-EFFECT structure variant
## (Stage-2b; independence family `random_walk / glmm-robust-re`)
##
## glmm_marginal and nimble_glmm both assume Gaussian random intercepts, so on
## well-behaved data they agree and the stacking weights collapse (they
## corroborate). This lens fits a genuinely DIFFERENT random-effects structure:
## the random intercepts are Student-t, b_g ~ t_nu(0, tau), so groups can take
## heavy-tailed (outlying) deviations. On data with an outlying group its
## predictive differs materially from the Gaussian-RE lenses', so LOO-stacking
## assigns it real weight --- the within-battery analogue of the population
## robust-t lens, and where koine's RE opinion is sensitive to random-effect
## structure rather than a weighted echo of one assumption.
##
## The t random effect is integrated out via its Normal scale-mixture
## representation (b | lambda ~ N(0, tau^2/lambda), lambda ~ Gamma(nu/2,nu/2)),
## so each Gauss-Laguerre node is an analytic Gaussian marginal (re_marginal.R).
## Sampled by DRAM over (beta, log sigma, log tau); it computes its OWN
## t-RE integrated predictive (stored natively), so assemble_loglik stacks the
## heavy-tailed-RE predictive rather than a Gaussian one (SS6.1). Gaussian
## response only (the analytic-per-node shortcut); nu is fixed.
## ----------------------------------------------------------------------------

# Unconstrained layout for the Gaussian robust-RE lens: coefficients, log-sigma,
# log-tau.
.glmm_robust_re_layout <- function(K) {
  list(dim = K + 2L, coef_idx = seq_len(K), sigma_idx = K + 1L,
       tau_idx = K + 2L)
}

# Log-posterior: the t-RE marginal log-likelihood + the Gaussian-family priors
# (shared with glmm_marginal, so its only difference is the random-effect
# distribution --- the diversity is structural, not a prior artefact).
.glmm_robust_re_logpost <- function(theta, design, redesign, prior, nu, mix,
                                    lay) {
  beta <- theta[lay$coef_idx]
  eta <- as.numeric(design$X %*% beta)
  log_sigma <- theta[lay$sigma_idx]
  log_tau <- theta[lay$tau_idx]
  ll <- .re_group_marginal_t(design$y, eta, redesign$group, exp(log_tau),
                             exp(log_sigma), nu, mix)
  ll + sum(stats::dnorm(beta, 0, prior$beta_sd, log = TRUE)) +
    stats::dnorm(log_sigma, 0, prior$log_sigma_sd, log = TRUE) +
    stats::dnorm(log_tau, 0, prior$log_tau_sd, log = TRUE)
}

#' Fit the heavy-tailed random-effect lens (Student-t random intercepts)
#'
#' Fits a Gaussian linear mixed model whose random intercepts are Student-t
#' (`b_g ~ t_nu(0, tau)`) by delayed-rejection adaptive Metropolis on the t-RE
#' marginal, and returns an [fb_fit()] carrying its own t-RE integrated
#' predictive. Resolves only for a Gaussian (identity) response with a supported
#' single random intercept; otherwise returns a lens-scope [abstention()].
#'
#' @param pass A `koine_pass` with a Gaussian family and a single `(1 | g)`.
#' @param lens_spec The lens's registry record; defaults to the `glmm_robust_re`
#'   record.
#' @param nu Degrees of freedom of the Student-t random effect (fixed; lower =
#'   heavier-tailed group effects). Default 4.
#' @param n_iter,n_burn,n_chains Sampler controls.
#' @param n_nodes Gauss-Laguerre nodes for the scale-mixture integral.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' g <- factor(rep(1:12, each = 8))
#' bt <- stats::rt(12, df = 4) * 1.2
#' d <- data.frame(x = rnorm(96), g = g)
#' d$y <- 1 + 0.5 * d$x + bt[g] + rnorm(96, 0, 0.7)
#' pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
#' fit_glmm_robust_re(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
#' }
fit_glmm_robust_re <- function(pass, lens_spec = NULL, nu = 4,
                               n_iter = 4000L, n_burn = 2000L, n_chains = 4L,
                               n_nodes = 32L) {
  if (is.null(lens_spec)) {
    lens_spec <- .registry_lens("glmm_robust_re")
  }
  if (!(pass$family$name == "gaussian" && pass$family$link == "identity")) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "glmm_robust_re serves a Gaussian (identity) response only."))
  }
  if (!isTRUE(pass$has_re) || !isTRUE(pass$re_supported)) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "glmm_robust_re requires a single supported random intercept."))
  }
  design <- .pass_design(pass)
  redesign <- .pass_re_design(pass)
  lay <- .glmm_robust_re_layout(design$K)
  prior <- .re_prior("gaussian")
  mix <- .t_mix_rule(nu, n_nodes)
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)

  start <- tryCatch(stats::coef(stats::lm.fit(design$X, design$y)),
                    error = function(e) rep.int(0, design$K))
  start[!is.finite(start)] <- 0
  logpost <- function(theta) {
    .glmm_robust_re_logpost(theta, design, redesign, prior, nu, mix, lay)
  }
  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    init <- numeric(lay$dim)
    init[lay$coef_idx] <- start + stats::rnorm(design$K, 0, 0.5)
    init[lay$sigma_idx] <- stats::rnorm(1, 0, 0.5)
    init[lay$tau_idx] <- stats::rnorm(1, 0, 0.5)
    chains[[c]] <- .dram_one_chain(init, logpost, n_iter = n_iter,
                                   n_burn = n_burn,
                                   adapt_start = max(100L, n_burn %/% 2L))
  }

  draws <- .glmm_robust_re_to_draws(chains, design, lay)
  log_lik <- .glmm_robust_re_loglik_array(chains, design, lay, nu, mix)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.0.0.9000", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(koine = "0.0.0.9000"),
    trace = list(fit = sprintf("glmm_robust_re (nu=%g, %d GL nodes) completed",
                               nu, n_nodes)))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, log_likelihood = log_lik, status = "contributed")
}

# (chain, draw, var) draws with sigma / tau back-transformed from the log scale.
.glmm_robust_re_to_draws <- function(chains, design, lay) {
  vnames <- c(paste0("b_", colnames(design$X)), "sigma", "tau")
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) {
    m <- chains[[c]]
    m[, lay$sigma_idx] <- exp(m[, lay$sigma_idx])
    m[, lay$tau_idx] <- exp(m[, lay$tau_idx])
    arr[, c, ] <- m
  }
  posterior::as_draws_array(arr)
}

# Native t-RE integrated predictive as array[chain, draw, obs] (the heavy-tailed
# random-effect marginalised out per observation).
.glmm_robust_re_loglik_array <- function(chains, design, lay, nu, mix) {
  cc <- length(chains)
  dd <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(cc, dd, design$N))
  for (c in seq_len(cc)) {
    m <- chains[[c]]
    for (s in seq_len(dd)) {
      eta <- as.numeric(design$X %*% m[s, lay$coef_idx])
      arr[c, s, ] <- .re_integrated_loglik_t(design$y, eta,
                                             tau = exp(m[s, lay$tau_idx]),
                                             sigma = exp(m[s, lay$sigma_idx]),
                                             nu = nu, mix = mix)
    }
  }
  arr
}
