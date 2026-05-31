## ----------------------------------------------------------------------------
## DRAM lens --- delayed-rejection adaptive Metropolis (koine-authored)
##
## A self-contained, base-R random-walk sampler in the `random_walk / dram`
## independence family. Adaptive Metropolis (Haario et al. 2001) gives the
## proposal covariance from chain history; a single delayed-rejection stage
## (Haario et al. 2006; Mira 2001) recovers acceptances a too-bold first
## proposal would have lost, using the exact second-stage acceptance ratio.
##
## The lens samples in an unconstrained parameter space (Gaussian dispersion is
## sampled as log-sigma) so the random walk is unconstrained; reported draws are
## back-transformed. It does not emit a pointwise log-likelihood --- koine
## computes that downstream via the family spec (F3).
## ----------------------------------------------------------------------------

# Draw one multivariate normal vector with covariance `sigma` (via Cholesky).
.rmvn <- function(mu, sigma) {
  L <- chol(sigma)
  as.numeric(mu + crossprod(L, stats::rnorm(length(mu))))
}

# Unconstrained parameter layout for a pass: coefficients first, then (for
# gaussian) a trailing log-sigma. Returns helpers to split / assemble.
.dram_param_layout <- function(pass, K) {
  has_sigma <- identical(pass$family$name, "gaussian")
  list(
    has_sigma = has_sigma,
    dim = K + as.integer(has_sigma),
    coef_idx = seq_len(K),
    sigma_idx = if (has_sigma) K + 1L else integer(0)
  )
}

# Log-posterior in the unconstrained space: family log-likelihood + weakly
# informative priors (beta ~ N(0, 5); log-sigma ~ N(0, 2)).
.dram_logpost <- function(theta, design, pass, lay, beta_prior_sd = 5) {
  beta <- theta[lay$coef_idx]
  eta <- as.numeric(design$X %*% beta)
  aux <- list()
  lp_extra <- 0
  if (lay$has_sigma) {
    log_sigma <- theta[lay$sigma_idx]
    aux$sigma <- exp(log_sigma)
    lp_extra <- stats::dnorm(log_sigma, 0, 2, log = TRUE)
  }
  ll <- sum(pass$family$loglik_pointwise(design$y, eta, aux))
  prior_beta <- sum(stats::dnorm(beta, 0, beta_prior_sd, log = TRUE))
  ll + prior_beta + lp_extra
}

# One chain of delayed-rejection adaptive Metropolis over an arbitrary log
# posterior `logpost(theta)`. Returns a (n_keep x dim) matrix of post-burn-in
# draws. Generic over the target so structure-variant lenses (e.g. robust-t)
# reuse the identical sampler.
.dram_one_chain <- function(init, logpost, n_iter, n_burn,
                            adapt_start, dr_scale = 0.2) {
  d <- length(init)
  x <- init
  lpx <- logpost(x)
  sig <- diag(0.1^2 / d, d)          # initial proposal covariance
  sd_scale <- 2.38^2 / d             # Haario optimal scaling
  eps <- 1e-8
  keep <- matrix(NA_real_, nrow = n_iter - n_burn, ncol = d)
  # running mean / covariance for adaptation
  run_mean <- x
  run_cov <- matrix(0, d, d)
  n_seen <- 1L

  for (i in seq_len(n_iter)) {
    # -- stage 1 proposal --
    y1 <- .rmvn(x, sig)
    lpy1 <- logpost(y1)
    a1 <- exp(min(0, lpy1 - lpx))
    accepted <- FALSE
    if (stats::runif(1) < a1) {
      x <- y1
      lpx <- lpy1
      accepted <- TRUE
    } else {
      # -- stage 2: delayed rejection with a shrunk proposal --
      sig2 <- sig * dr_scale^2
      y2 <- .rmvn(x, sig2)
      lpy2 <- logpost(y2)
      # exact DR stage-2 acceptance (symmetric stage-1 kernel cancels in the
      # q-ratio except through the y1 conditional; see Haario et al. 2006).
      a1_y2 <- exp(min(0, lpy1 - lpy2))
      log_q_num <- .dmvn_log(y1, y2, sig)
      log_q_den <- .dmvn_log(y1, x, sig)
      num <- lpy2 + log_q_num + log1p(-a1_y2)
      den <- lpx + log_q_den + log1p(-a1)
      a2 <- if (a1 >= 1 || a1_y2 >= 1) 0 else exp(min(0, num - den))
      if (is.finite(a2) && stats::runif(1) < a2) {
        x <- y2
        lpx <- lpy2
        accepted <- TRUE
      }
    }
    # -- adaptation (recursive mean/cov; engaged after adapt_start) --
    n_seen <- n_seen + 1L
    prev_mean <- run_mean
    run_mean <- prev_mean + (x - prev_mean) / n_seen
    run_cov <- run_cov +
      (tcrossprod(x - prev_mean, x - run_mean) - run_cov) / n_seen
    if (i >= adapt_start) {
      sig <- sd_scale * (run_cov + eps * diag(d))
    }
    if (i > n_burn) keep[i - n_burn, ] <- x
  }
  keep
}

# Log density of a multivariate normal N(mu, sigma) at x (for the DR q-ratio).
.dmvn_log <- function(x, mu, sigma) {
  d <- length(x)
  L <- chol(sigma)
  z <- backsolve(L, x - mu, transpose = TRUE)
  -0.5 * (d * log(2 * pi)) - sum(log(diag(L))) - 0.5 * sum(z^2)
}

#' Fit the DRAM lens
#'
#' Runs delayed-rejection adaptive Metropolis on the pass's GLM and returns an
#' [fb_fit()]. The fit carries equal-weight posterior draws and a complete
#' provenance record; the pointwise log-likelihood is left `NULL` for
#' [assemble_loglik()] to compute via the family spec (F3).
#'
#' @param pass A `koine_pass`.
#' @param lens_spec The lens's registry record (for provenance + the derived
#'   seed); defaults to the registry's `dram` record.
#' @param n_iter Iterations per chain.
#' @param n_burn Burn-in iterations discarded per chain.
#' @param n_chains Number of chains.
#' @return A `koine_fb_fit`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(40))
#' d$y <- 1 + 0.8 * d$x + rnorm(40)
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
#' }
#' @export
fit_dram <- function(pass, lens_spec = NULL, n_iter = 4000L, n_burn = 2000L,
                     n_chains = 4L) {
  if (is.null(lens_spec)) {
    reg <- koine_stage1_registry()
    lens_spec <- reg[[which(vapply(reg, function(z) z$id, "") == "dram")]]
  }
  design <- .pass_design(pass)
  lay <- .dram_param_layout(pass, design$K)
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)

  # dispersed inits: OLS/GLM-ish start jittered per chain.
  # Warm-start init only; its convergence notices on awkward data are noise
  # (the MCMC sampler is the inference), so they are suppressed, not surfaced.
  start <- tryCatch(
    suppressWarnings(stats::coef(stats::glm.fit(
      design$X, design$y, family = .pass_glm_family(pass)))),
    error = function(e) rep.int(0, design$K)
  )
  start[!is.finite(start)] <- 0
  logpost <- function(theta) .dram_logpost(theta, design, pass, lay)
  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    init <- numeric(lay$dim)
    init[lay$coef_idx] <- start + stats::rnorm(design$K, 0, 0.5)
    if (lay$has_sigma) init[lay$sigma_idx] <- stats::rnorm(1, 0, 0.5)
    chains[[c]] <- .dram_one_chain(init, logpost,
                                   n_iter = n_iter, n_burn = n_burn,
                                   adapt_start = max(100L, n_burn %/% 2L))
  }

  draws <- .dram_to_draws(chains, design, lay)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.0.0.9000", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(koine = "0.0.0.9000"),
    trace = list(fit = "dram completed")
  )
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

# Map a base-R family onto a glm.fit family for the warm-start only.
.pass_glm_family <- function(pass) {
  switch(pass$family$name,
         gaussian = stats::gaussian(link = pass$family$link),
         binomial = stats::binomial(link = pass$family$link),
         poisson = stats::poisson(link = pass$family$link),
         stats::gaussian())
}

# Assemble chains into a posterior draws_array with named variables
# (b_<col> for coefficients, sigma for the gaussian dispersion).
.dram_to_draws <- function(chains, design, lay) {
  vnames <- paste0("b_", colnames(design$X))
  if (lay$has_sigma) vnames <- c(vnames, "sigma")
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) {
    m <- chains[[c]]
    if (lay$has_sigma) m[, lay$sigma_idx] <- exp(m[, lay$sigma_idx])
    arr[, c, ] <- m
  }
  posterior::as_draws_array(arr)
}
