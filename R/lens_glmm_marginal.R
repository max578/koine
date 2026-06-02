## ----------------------------------------------------------------------------
## glmm_marginal lens --- Hamiltonian (NUTS) MCMC on the marginal mixed model
## (Stage-2; independence family `hamiltonian / glmm-marginal-quadrature`)
##
## A koine-authored mixed-model lens that *never samples the random effects*: it
## integrates them out and runs the koine-authored No-U-Turn Sampler
## (`.nuts_one_chain`) over `(beta, log sigma, log tau)`, the fixed coefficients
## plus the residual and random-intercept SDs. The per-group marginal likelihood
## is exact for the Gaussian LMM (a multivariate normal per group) and adaptive
## Gauss-Hermite quadrature for a GLMM (`.re_group_marginal_loglik`); its
## gradient --- closed form for the Gaussian, the score-identity quadrature
## estimate for a GLMM (`.re_group_marginal_grad`) --- drives the Hamiltonian
## trajectory.
##
## The kernel is gradient-informed by design: the Poisson marginal's confounded
## count-intercept ridge mixes far too slowly under a random-walk kernel (an
## SBC-failing intercept even after the level reparameterisation), but a single
## Hamiltonian trajectory traverses it. This is deliberately orthogonal to the
## sealed engines: NIMBLE *samples* the random effects (Gibbs), INLA
## *Laplace-approximates* the integral, and this lens *quadrature-integrates* it
## under a Hamiltonian kernel --- three distinct treatments of the one integral,
## so their agreement (or disagreement) carries a genuine fourth-witness signal,
## strongest for GLMMs where the Laplace approximation is only approximate. It
## reports `(beta, sigma, tau)` draws and leaves the pointwise log-likelihood for
## `assemble_loglik` to compute as the per-observation integrated predictive
## (F3), so it is stacked on the identical held-out target as the NIMBLE GLMM
## lens.
## ----------------------------------------------------------------------------

# Unconstrained parameter layout: coefficients, then (Gaussian) a log-sigma,
# then always a trailing log-tau (the random-intercept SD, sampled on the log
# scale so the Hamiltonian trajectory is unconstrained and tau stays positive).
.glmm_marginal_layout <- function(pass, K) {
  has_sigma <- identical(pass$family$name, "gaussian")
  list(
    has_sigma = has_sigma,
    dim = K + as.integer(has_sigma) + 1L,
    coef_idx = seq_len(K),
    sigma_idx = if (has_sigma) K + 1L else integer(0),
    tau_idx = K + as.integer(has_sigma) + 1L
  )
}

# Value AND gradient of the marginal log-posterior in the unconstrained space
# (theta -> list(value, grad)), the differentiable target the NUTS kernel rides.
# The gradient assembles three pieces from `.re_group_marginal_grad()`:
#   * coefficients: d/d beta = X' s_obs + Normal-prior score (s_obs is the
#     expected per-observation score, exact for the Gaussian LMM and the
#     score-identity quadrature estimate for a GLMM);
#   * log sigma (Gaussian): the marginal's log-sigma derivative + prior;
#   * log tau: the random-effect-density derivative + prior, plus --- under the
#     log-link Poisson level reparameterisation (`rep_pos`) --- the chain term
#     from beta0 = alpha0 - tau^2 / 2, which couples the intercept coordinate's
#     full (data + prior) gradient into log tau by exactly -tau^2. Verified
#     against finite differences in `test-glmm_marginal`.
.glmm_marginal_logpost_grad <- function(theta, design, redesign, pass, lay, gh,
                                        prior, rep_pos = 0L) {
  beta <- theta[lay$coef_idx]
  log_tau <- theta[lay$tau_idx]
  tau <- exp(log_tau)
  if (rep_pos > 0L) {
    beta[rep_pos] <- beta[rep_pos] - tau^2 / 2
  }
  eta_fixed <- as.numeric(design$X %*% beta)
  aux <- list()
  lp_extra <- stats::dnorm(log_tau, 0, prior$log_tau_sd, log = TRUE)
  if (lay$has_sigma) {
    log_sigma <- theta[lay$sigma_idx]
    aux$sigma <- exp(log_sigma)
    lp_extra <- lp_extra + stats::dnorm(log_sigma, 0, prior$log_sigma_sd,
                                        log = TRUE)
  }
  mg <- .re_group_marginal_grad(pass$family, design$y, eta_fixed,
                                redesign$group, tau, aux, gh)
  prior_beta <- sum(stats::dnorm(beta, 0, prior$beta_sd, log = TRUE))
  value <- mg$value + prior_beta + lp_extra
  if (!is.finite(value)) {
    return(list(value = -Inf, grad = numeric(lay$dim)))
  }
  grad <- numeric(lay$dim)
  xts <- as.numeric(crossprod(design$X, mg$s_obs))
  prior_score_beta <- -beta / prior$beta_sd^2
  grad[lay$coef_idx] <- xts + prior_score_beta
  grad[lay$tau_idx] <- mg$d_logtau - log_tau / prior$log_tau_sd^2
  if (rep_pos > 0L) {
    # beta0 = alpha0 - tau^2/2 : d/d log tau picks up -tau^2 times the full
    # intercept-coordinate gradient (data score + prior score).
    grad[lay$tau_idx] <- grad[lay$tau_idx] -
      tau^2 * (xts[rep_pos] + prior_score_beta[rep_pos])
  }
  if (lay$has_sigma) {
    grad[lay$sigma_idx] <- mg$d_logsigma -
      theta[lay$sigma_idx] / prior$log_sigma_sd^2
  }
  if (!all(is.finite(grad))) {
    return(list(value = -Inf, grad = numeric(lay$dim)))
  }
  list(value = value, grad = grad)
}

#' Fit the marginal mixed-model lens
#'
#' Runs the koine-authored No-U-Turn Sampler on the random-effects *marginal*
#' model --- the random intercepts integrated out, exactly for the Gaussian LMM
#' and by adaptive Gauss-Hermite quadrature for a GLMM --- and returns an
#' [fb_fit()] of `(beta, sigma, tau)` draws. The gradient-informed kernel mixes
#' the confounded count-intercept ridge of the Poisson marginal that a
#' random-walk kernel could not. The pointwise log-likelihood is left `NULL` for
#' [assemble_loglik()] to compute as the per-observation integrated predictive
#' (F3). If the pass carries no supported random intercept, returns a lens-scope
#' [abstention()] rather than a wrong fit (I1).
#'
#' @param pass A `koine_pass` carrying a single random intercept `(1 | g)`.
#' @param lens_spec The lens's registry record; defaults to the registry's
#'   `glmm_marginal` record.
#' @param n_iter Iterations per chain (including warm-up). The Hamiltonian kernel
#'   needs far fewer than a random walk for the same effective sample size, so
#'   the default is modest.
#' @param n_burn Warm-up iterations per chain (step-size + mass-matrix
#'   adaptation), discarded from the returned draws.
#' @param n_chains Number of chains.
#' @param n_nodes Gauss-Hermite nodes for the GLMM marginal (ignored for the
#'   exact Gaussian LMM path).
#' @param max_treedepth Maximum NUTS tree depth, capping the leapfrog budget per
#'   iteration at `2^max_treedepth`. Bounds the cost on a hard, wide-tau marginal
#'   (whose long trajectories would otherwise run away); saturation of the cap is
#'   surfaced in the provenance trace, never silently absorbed.
#' @return A `koine_fb_fit`, or a `koine_abstention` for an unsupported pass.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' g <- factor(rep(1:8, each = 5))
#' d <- data.frame(x = rnorm(40), g = g)
#' d$y <- 1 + 0.5 * d$x + rnorm(8)[g] + rnorm(40, 0, 0.7)
#' pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
#' fit_glmm_marginal(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
#' }
fit_glmm_marginal <- function(pass, lens_spec = NULL, n_iter = 1500L,
                              n_burn = 750L, n_chains = 4L, n_nodes = 20L,
                              max_treedepth = 9L) {
  if (is.null(lens_spec)) {
    lens_spec <- .registry_lens("glmm_marginal")
  }
  if (!isTRUE(pass$has_re) || !isTRUE(pass$re_supported)) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "glmm_marginal requires a single supported random intercept."))
  }
  design <- .pass_design(pass)
  redesign <- .pass_re_design(pass)
  lay <- .glmm_marginal_layout(pass, design$K)
  gh <- .gh_rule(n_nodes)
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)

  start <- tryCatch(
    suppressWarnings(stats::coef(stats::glm.fit(
      design$X, design$y, family = .pass_glm_family(pass)))),
    error = function(e) rep.int(0, design$K))
  start[!is.finite(start)] <- 0
  prior <- .re_prior(pass$family$name)
  # log-link Poisson: sample the decorrelated marginal-mean-level intercept (the
  # population GLM intercept already estimates it, so the warm start is apt).
  rep_pos <- if (identical(pass$family$name, "poisson") &&
                   identical(pass$family$link, "log")) {
    match("(Intercept)", colnames(design$X), nomatch = 0L)
  } else {
    0L
  }
  lpg <- function(theta) {
    .glmm_marginal_logpost_grad(theta, design, redesign, pass, lay, gh, prior,
                                rep_pos)
  }
  chains <- vector("list", n_chains)
  divergences <- 0L
  saturated <- 0L
  for (c in seq_len(n_chains)) {
    init <- numeric(lay$dim)
    init[lay$coef_idx] <- start + stats::rnorm(design$K, 0, 0.5)
    if (lay$has_sigma) init[lay$sigma_idx] <- stats::rnorm(1, 0, 0.5)
    init[lay$tau_idx] <- stats::rnorm(1, 0, 0.5)
    # target_accept 0.9 (above the NUTS default 0.8): the weakly-identified
    # variance component gives the marginal a tau funnel, and the smaller adapted
    # step size clears the neck with fewer divergences. max_treedepth caps the
    # leapfrog budget per iteration (2^depth) so a hard, wide-tau marginal cannot
    # run unbounded; saturation of the cap is surfaced (never silently absorbed).
    ch <- .nuts_one_chain(init, lpg, n_iter = n_iter, n_warmup = n_burn,
                          target_accept = 0.9, max_treedepth = max_treedepth)
    chains[[c]] <- ch$draws
    divergences <- divergences + ch$divergences
    saturated <- saturated + ch$saturated
  }

  draws <- .glmm_marginal_to_draws(chains, design, lay, rep_pos)
  # The Gaussian LMM marginal is exact; the GLMM marginal is a quadrature
  # approximation --- reported honestly per family (I1).
  exact <- if (lay$has_sigma) "exact" else "approximate"
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.0.0.9000", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = exact,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(koine = "0.0.0.9000"),
    trace = list(fit = sprintf(
      "glmm_marginal/NUTS completed (%d GH nodes, %d divergences, %d treedepth-saturated)",
      n_nodes, divergences, saturated)))
  fb_fit(draws = draws, provenance = prov, exactness = exact,
         sample_stats = NULL, status = "contributed")
}

# Assemble chains into a draws_array with named variables (b_<col> for
# coefficients, sigma for the Gaussian residual SD, tau for the random-intercept
# SD), back-transforming the log-scale sigma / tau and (when `rep_pos` > 0)
# recovering the coefficient intercept beta0 = alpha0 - tau^2 / 2 from its
# sampled marginal-mean-level form.
.glmm_marginal_to_draws <- function(chains, design, lay, rep_pos = 0L) {
  vnames <- paste0("b_", colnames(design$X))
  if (lay$has_sigma) vnames <- c(vnames, "sigma")
  vnames <- c(vnames, "tau")
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) {
    m <- chains[[c]]
    if (lay$has_sigma) m[, lay$sigma_idx] <- exp(m[, lay$sigma_idx])
    tau <- exp(m[, lay$tau_idx])
    m[, lay$tau_idx] <- tau
    if (rep_pos > 0L) {
      m[, rep_pos] <- m[, rep_pos] - tau^2 / 2
    }
    arr[, c, ] <- m
  }
  posterior::as_draws_array(arr)
}
