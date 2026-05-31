## ----------------------------------------------------------------------------
## glmm_marginal lens --- random-walk MCMC on the marginal mixed model
## (Stage-2; independence family `random_walk / glmm-marginal-quadrature`)
##
## A koine-authored mixed-model lens that *never samples the random effects*: it
## integrates them out and runs the existing delayed-rejection adaptive
## Metropolis sampler (`.dram_one_chain`) over `(beta, log sigma, log tau)`, the
## fixed coefficients plus the residual and random-intercept SDs.
## The per-group marginal likelihood is exact for the Gaussian LMM (a
## multivariate normal per group) and Gauss-Hermite quadrature for a GLMM
## (`.re_group_marginal_loglik`).
##
## This is deliberately orthogonal to the sealed engines: NIMBLE *samples* the
## random effects (Gibbs), INLA *Laplace-approximates* the integral, and this
## lens *quadrature-integrates* it under a random-walk kernel --- three distinct
## treatments of the one integral, so their agreement (or disagreement) carries
## a genuine fourth-witness signal, strongest for GLMMs where the Laplace
## approximation is only approximate. It reports `(beta, sigma, tau)` draws and
## leaves the pointwise log-likelihood for `assemble_loglik` to compute as the
## per-observation integrated predictive (F3), so it is stacked on the identical
## held-out target as the NIMBLE GLMM lens.
## ----------------------------------------------------------------------------

# Unconstrained parameter layout: coefficients, then (Gaussian) a log-sigma,
# then always a trailing log-tau (the random-intercept SD, sampled on the log
# scale so the random walk is unconstrained and tau stays positive).
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

# Log-posterior in the unconstrained space: the random-effects marginal
# log-likelihood + weakly-informative priors (family-appropriate SDs from
# `.re_prior()`: coefficients, the Gaussian residual log-sigma, and the
# random-intercept log-tau). The priors match the NIMBLE lens so cross-lens
# agreement is method-independent rather than a prior artefact.
.glmm_marginal_logpost <- function(theta, design, redesign, pass, lay, gh,
                                   prior) {
  beta <- theta[lay$coef_idx]
  eta_fixed <- as.numeric(design$X %*% beta)
  log_tau <- theta[lay$tau_idx]
  tau <- exp(log_tau)
  aux <- list()
  lp_extra <- stats::dnorm(log_tau, 0, prior$log_tau_sd, log = TRUE)
  if (lay$has_sigma) {
    log_sigma <- theta[lay$sigma_idx]
    aux$sigma <- exp(log_sigma)
    lp_extra <- lp_extra + stats::dnorm(log_sigma, 0, prior$log_sigma_sd,
                                        log = TRUE)
  }
  ll <- .re_group_marginal_loglik(pass$family, design$y, eta_fixed,
                                  redesign$group, tau, aux, gh)
  prior_beta <- sum(stats::dnorm(beta, 0, prior$beta_sd, log = TRUE))
  ll + prior_beta + lp_extra
}

#' Fit the marginal mixed-model lens
#'
#' Runs delayed-rejection adaptive Metropolis on the random-effects *marginal*
#' model --- the random intercepts integrated out, exactly for the Gaussian LMM
#' and by Gauss-Hermite quadrature for a GLMM --- and returns an [fb_fit()] of
#' `(beta, sigma, tau)` draws. The pointwise log-likelihood is left `NULL` for
#' [assemble_loglik()] to compute as the per-observation integrated predictive
#' (F3). If the pass carries no supported random intercept, returns a lens-scope
#' [abstention()] rather than a wrong fit (I1).
#'
#' @param pass A `koine_pass` carrying a single random intercept `(1 | g)`.
#' @param lens_spec The lens's registry record; defaults to the registry's
#'   `glmm_marginal` record.
#' @param n_iter Iterations per chain (including burn-in).
#' @param n_burn Burn-in iterations discarded per chain.
#' @param n_chains Number of chains.
#' @param n_nodes Gauss-Hermite nodes for the GLMM marginal (ignored for the
#'   exact Gaussian LMM path).
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
fit_glmm_marginal <- function(pass, lens_spec = NULL, n_iter = 4000L,
                              n_burn = 2000L, n_chains = 4L, n_nodes = 20L) {
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
  logpost <- function(theta) {
    .glmm_marginal_logpost(theta, design, redesign, pass, lay, gh, prior)
  }
  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    init <- numeric(lay$dim)
    init[lay$coef_idx] <- start + stats::rnorm(design$K, 0, 0.5)
    if (lay$has_sigma) init[lay$sigma_idx] <- stats::rnorm(1, 0, 0.5)
    init[lay$tau_idx] <- stats::rnorm(1, 0, 0.5)
    chains[[c]] <- .dram_one_chain(init, logpost, n_iter = n_iter,
                                   n_burn = n_burn,
                                   adapt_start = max(100L, n_burn %/% 2L))
  }

  draws <- .glmm_marginal_to_draws(chains, design, lay)
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
    trace = list(fit = sprintf("glmm_marginal completed (%d GH nodes)",
                               n_nodes)))
  fb_fit(draws = draws, provenance = prov, exactness = exact,
         sample_stats = NULL, status = "contributed")
}

# Assemble chains into a draws_array with named variables (b_<col> for
# coefficients, sigma for the Gaussian residual SD, tau for the random-intercept
# SD), back-transforming the log-scale sigma / tau.
.glmm_marginal_to_draws <- function(chains, design, lay) {
  vnames <- paste0("b_", colnames(design$X))
  if (lay$has_sigma) vnames <- c(vnames, "sigma")
  vnames <- c(vnames, "tau")
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) {
    m <- chains[[c]]
    if (lay$has_sigma) m[, lay$sigma_idx] <- exp(m[, lay$sigma_idx])
    m[, lay$tau_idx] <- exp(m[, lay$tau_idx])
    arr[, c, ] <- m
  }
  posterior::as_draws_array(arr)
}
