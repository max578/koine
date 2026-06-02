## ----------------------------------------------------------------------------
## NIMBLE GLMM lens (Stage-2; independence family `gibbs / nimble-mwg-glmm`)
##
## NIMBLE's conjugate / Metropolis-within-Gibbs MCMC on a generalised linear
## mixed model with a single random intercept: it *samples* the random effects
## b[g] ~ N(0, tau^2) explicitly. This is the orthogonal counterpart to the
## `glmm_marginal` lens, which integrates the same effects out --- the two
## carry identical weakly-informative priors (beta ~ N(0, 5); log-sigma,
## log-tau ~ N(0, 2)), so any agreement between them is method-independent
## consensus and any disagreement is a genuine inference-family signal (SS2).
## Sampling b (Gibbs) is also orthogonal to INLA's Laplace marginalisation, the
## sealed engine koine must not echo.
##
## The lens monitors only `(beta, sigma, tau)` --- it discards the sampled
## effects, because the fourth opinion is stacked on the per-observation
## integrated predictive (random effect marginalised under the sampled tau,
## `assemble_loglik`), the identical held-out target the marginal lens uses.
## ----------------------------------------------------------------------------

# Build the family-specific NIMBLE GLMM code (single random intercept), or NULL
# if the family/link pair is unsupported (the caller then abstains, lens scope).
.nimble_glmm_code <- function(family_name, link) {
  ok <- (family_name == "gaussian" && link == "identity") ||
    (family_name == "binomial" && link == "logit") ||
    (family_name == "poisson" && link == "log")
  if (!ok) return(NULL)
  # Prior SDs (beta_sd, log_sigma_sd, log_tau_sd) are supplied as model
  # constants by the caller (`.re_prior()`), so the family-appropriate prior is
  # applied without duplicating the model code per prior.
  switch(
    family_name,
    gaussian = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = beta_sd)
      log_sigma ~ dnorm(0, sd = log_sigma_sd)
      sigma <- exp(log_sigma)
      log_tau ~ dnorm(0, sd = log_tau_sd)
      tau <- exp(log_tau)
      for (j in 1:G) b[j] ~ dnorm(0, sd = tau)
      for (i in 1:N) {
        mu[i] <- inprod(X[i, 1:K], beta[1:K]) + b[group[i]]
        y[i] ~ dnorm(mu[i], sd = sigma)
      }
    }),
    binomial = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = beta_sd)
      log_tau ~ dnorm(0, sd = log_tau_sd)
      tau <- exp(log_tau)
      for (j in 1:G) b[j] ~ dnorm(0, sd = tau)
      for (i in 1:N) {
        logit(p[i]) <- inprod(X[i, 1:K], beta[1:K]) + b[group[i]]
        y[i] ~ dbin(p[i], 1)
      }
    }),
    poisson = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = beta_sd)
      log_tau ~ dnorm(0, sd = log_tau_sd)
      tau <- exp(log_tau)
      for (j in 1:G) b[j] ~ dnorm(0, sd = tau)
      for (i in 1:N) {
        log(lambda[i]) <- inprod(X[i, 1:K], beta[1:K]) + b[group[i]]
        y[i] ~ dpois(lambda[i])
      }
    })
  )
}

#' Fit the NIMBLE GLMM lens
#'
#' Runs NIMBLE's conjugate / Metropolis-within-Gibbs sampler on a mixed model
#' with a single random intercept and returns an [fb_fit()] of `(beta, sigma,
#' tau)` draws --- the random effects are sampled but not retained, since the
#' opinion is stacked on the integrated predictive. If the family/link pair is
#' unsupported, or the pass carries no supported random intercept, returns a
#' lens-scope [abstention()] rather than a wrong fit (I1).
#'
#' @param pass A `koine_pass` carrying a single random intercept `(1 | g)`.
#' @param lens_spec The lens's registry record; defaults to the registry's
#'   `nimble_glmm` record.
#' @param n_iter Iterations per chain (including burn-in).
#' @param n_burn Burn-in iterations discarded per chain.
#' @param n_chains Number of chains.
#' @return A `koine_fb_fit`, or a `koine_abstention` for an unsupported pass.
#' @examplesIf requireNamespace("nimble", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' g <- factor(rep(1:8, each = 5))
#' d <- data.frame(x = rnorm(40), g = g)
#' d$y <- 1 + 0.5 * d$x + rnorm(8)[g] + rnorm(40, 0, 0.7)
#' pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
#' fit_nimble_glmm(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_nimble_glmm <- function(pass, lens_spec = NULL, n_iter = 4000L,
                            n_burn = 2000L, n_chains = 4L) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the nimble package is not installed."))
  }
  if (!"package:nimble" %in% search()) {
    suppressMessages(attachNamespace("nimble"))
  }
  if (is.null(lens_spec)) {
    lens_spec <- .registry_lens("nimble_glmm")
  }
  if (!isTRUE(pass$has_re) || !isTRUE(pass$re_supported)) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "nimble_glmm requires a single supported random intercept."))
  }
  code <- .nimble_glmm_code(pass$family$name, pass$family$link)
  if (is.null(code)) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("nimble_glmm has no model for family '%s' + link '%s'.",
                       pass$family$name, pass$family$link)))
  }

  design <- .pass_design(pass)
  redesign <- .pass_re_design(pass)
  has_sigma <- identical(pass$family$name, "gaussian")
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  monitors <- c("beta", "tau")
  if (has_sigma) monitors <- c(monitors, "sigma")

  prior <- .re_prior(pass$family$name)
  constants <- list(N = design$N, K = design$K, X = design$X,
                    G = redesign$n_groups, group = redesign$group,
                    beta_sd = prior$beta_sd, log_tau_sd = prior$log_tau_sd)
  if (has_sigma) constants$log_sigma_sd <- prior$log_sigma_sd
  data_list <- list(y = design$y)
  warm <- tryCatch(
    suppressWarnings(stats::coef(stats::glm.fit(
      design$X, design$y, family = .pass_glm_family(pass)))),
    error = function(e) rep.int(0, design$K))
  warm[!is.finite(warm)] <- 0
  inits <- list(beta = as.numeric(warm), log_tau = 0,
                b = rep.int(0, redesign$n_groups))
  if (has_sigma) inits$log_sigma <- 0

  samples <- NULL
  .nimble_quiet({
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = data_list, inits = inits,
                                 calculate = FALSE)
    cmodel <- nimble::compileNimble(model)
    conf <- nimble::configureMCMC(model, monitors = monitors, print = FALSE)
    mcmc <- nimble::buildMCMC(conf)
    cmcmc <- nimble::compileNimble(mcmc, project = model)
    samples <- nimble::runMCMC(
      cmcmc, niter = n_iter, nburnin = n_burn, nchains = n_chains,
      setSeed = dseed + seq_len(n_chains), progressBar = FALSE,
      samplesAsCodaMCMC = FALSE)
  })
  if (!is.list(samples)) samples <- list(samples)

  draws <- .nimble_glmm_to_draws(samples, design, has_sigma)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.1.0", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(nimble = as.character(utils::packageVersion("nimble"))),
    trace = list(fit = "nimble GLMM MCMC completed"))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

# Validation-only NIMBLE reference for the heavy-tailed random-effect lens: a
# Gaussian LMM with Student-t random intercepts written as the SAME Normal
# scale-mixture the marginal lens uses (b_j ~ N(0, tau^2 / lam_j), lam_j ~
# Gamma(nu/2, nu/2)), so the two are concordance-comparable on the identical
# t-RE definition. Not a registry lens --- it is the independent reference the
# koine-authored glmm_robust_re lens is validated against (SS11 concordance).
.fit_nimble_glmm_t <- function(pass, nu = 4, n_iter = 4000L, n_burn = 2000L,
                               n_chains = 3L) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the nimble package is not installed."))
  }
  if (!"package:nimble" %in% search()) {
    suppressMessages(attachNamespace("nimble"))
  }
  design <- .pass_design(pass)
  redesign <- .pass_re_design(pass)
  dseed <- .derived_seed(pass$seed, "nimble_glmm_t")
  prior <- .re_prior("gaussian")
  code <- nimble::nimbleCode({
    for (k in 1:K) beta[k] ~ dnorm(0, sd = beta_sd)
    log_sigma ~ dnorm(0, sd = log_sigma_sd)
    sigma <- exp(log_sigma)
    log_tau ~ dnorm(0, sd = log_tau_sd)
    tau <- exp(log_tau)
    for (j in 1:G) {
      lam[j] ~ dgamma(nu_half, nu_half)
      b[j] ~ dnorm(0, sd = tau / sqrt(lam[j]))
    }
    for (i in 1:N) {
      mu[i] <- inprod(X[i, 1:K], beta[1:K]) + b[group[i]]
      y[i] ~ dnorm(mu[i], sd = sigma)
    }
  })
  constants <- list(N = design$N, K = design$K, X = design$X,
                    G = redesign$n_groups, group = redesign$group,
                    beta_sd = prior$beta_sd, log_sigma_sd = prior$log_sigma_sd,
                    log_tau_sd = prior$log_tau_sd, nu_half = nu / 2)
  warm <- tryCatch(stats::coef(stats::lm.fit(design$X, design$y)),
                   error = function(e) rep.int(0, design$K))
  warm[!is.finite(warm)] <- 0
  inits <- list(beta = as.numeric(warm), log_sigma = 0, log_tau = 0,
                lam = rep.int(1, redesign$n_groups),
                b = rep.int(0, redesign$n_groups))
  samples <- NULL
  .nimble_quiet({
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = list(y = design$y), inits = inits,
                                 calculate = FALSE)
    cmodel <- nimble::compileNimble(model)
    conf <- nimble::configureMCMC(model, monitors = c("beta", "sigma", "tau"),
                                  print = FALSE)
    cmcmc <- nimble::compileNimble(nimble::buildMCMC(conf), project = model)
    samples <- nimble::runMCMC(cmcmc, niter = n_iter, nburnin = n_burn,
                               nchains = n_chains,
                               setSeed = dseed + seq_len(n_chains),
                               progressBar = FALSE, samplesAsCodaMCMC = FALSE)
  })
  if (!is.list(samples)) samples <- list(samples)
  prov <- provenance_record(
    method = "nimble_glmm_t", method_version = "0.1.0",
    impl_language = "R", seed = pass$seed, derived_seed = dseed,
    ir_hash = pass$ir_hash, maturity_tier = "reference-validated",
    independence_key = list(family = "gibbs", kernel = "nimble-mwg-glmm-t"),
    role = "opinion", exactness = "exact",
    lib_versions = list(nimble = as.character(utils::packageVersion("nimble"))),
    trace = list(fit = "nimble t-RE reference completed"))
  draws <- .nimble_glmm_to_draws(samples, design, TRUE)
  fb_fit(draws = draws, provenance = prov, exactness = "exact",
         sample_stats = NULL, status = "contributed")
}

# Convert NIMBLE's per-chain sample matrices (columns beta[1].., tau, [sigma])
# into a draws_array with the b_<col> / sigma / tau naming the assembly expects.
.nimble_glmm_to_draws <- function(samples, design, has_sigma) {
  bnames <- paste0("b_", colnames(design$X))
  target <- c(bnames, if (has_sigma) "sigma", "tau")
  n_keep <- nrow(samples[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(samples), length(target)),
               dimnames = list(NULL, NULL, target))
  for (c in seq_along(samples)) {
    m <- samples[[c]]
    src <- matrix(NA_real_, nrow = n_keep, ncol = length(target),
                  dimnames = list(NULL, target))
    for (k in seq_len(design$K)) {
      src[, bnames[k]] <- m[, paste0("beta[", k, "]")]
    }
    if (has_sigma) src[, "sigma"] <- m[, "sigma"]
    src[, "tau"] <- m[, "tau"]
    arr[, c, ] <- src
  }
  posterior::as_draws_array(arr)
}
