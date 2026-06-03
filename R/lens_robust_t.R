## ----------------------------------------------------------------------------
## Robust-t lens --- a MODEL-STRUCTURE variant (the diversity engine)
##
## Unlike DRAM / NIMBLE / PG-Gibbs / DREAM, which are different *inference
## families* on the *same* posterior (and therefore corroborate, collapsing the
## stacking weights), the robust-t lens fits a genuinely DIFFERENT model: a
## Student-t likelihood (heavy tails) in place of the Gaussian. On data with
## outliers its predictive distribution differs materially from the Gaussian
## lenses', so LOO-stacking assigns it real weight --- this is where koine's
## orthogonal-diversity premise becomes visible: the disagreement carries a new
## signal ("sensitive to model structure"), not just a consistency check.
##
## It is sampled by the same delayed-rejection adaptive Metropolis kernel as
## DRAM, and it computes its OWN pointwise log-likelihood under the t model
## (stored natively in the fb_fit), so assemble_loglik stacks its t-predictive
## rather than a Gaussian one (SS6.1).
## ----------------------------------------------------------------------------

# Scaled Student-t log density: log [ (1/sigma) * t_nu((y - mu)/sigma) ].
.student_t_loglik <- function(y, mu, sigma, nu) {
  stats::dt((y - mu) / sigma, df = nu, log = TRUE) - log(sigma)
}

#' Fit the robust-t lens (Student-t linear model)
#'
#' Fits `y ~ Student-t(nu, X beta, sigma)` by delayed-rejection adaptive
#' Metropolis (the DRAM kernel) and returns an [fb_fit()] carrying its own
#' Student-t pointwise log-likelihood. Resolves only for a continuous
#' (Gaussian-declared) response; otherwise returns a lens-scope [abstention()].
#'
#' @param pass A `koine_pass` with a Gaussian (identity) family.
#' @param lens_spec The lens's registry record; defaults to the `robust_t`
#'   record.
#' @param nu Degrees of freedom of the Student-t likelihood (fixed; lower = of
#'   heavier tails). Default 4.
#' @param n_iter,n_burn,n_chains Sampler controls.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(50))
#' d$y <- 1 + 0.8 * d$x + stats::rt(50, df = 4) * 0.6
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit_robust_t(pass, nu = 4, n_iter = 800L, n_burn = 400L, n_chains = 2L)
#' }
#' @export
fit_robust_t <- function(pass, lens_spec = NULL, nu = 4,
                         n_iter = 4000L, n_burn = 2000L, n_chains = 4L) {
  if (!(pass$family$name == "gaussian" && pass$family$link == "identity")) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "robust_t serves a continuous (Gaussian) response only."))
  }
  if (is.null(lens_spec)) {
    reg <- koine_stage1_registry()
    lens_spec <- reg[[which(vapply(reg, function(z) z$id, "") == "robust_t")]]
  }
  design <- .pass_design(pass)
  K <- design$K
  coef_idx <- seq_len(K)
  sigma_idx <- K + 1L

  logpost <- function(theta) {
    beta <- theta[coef_idx]
    ls <- theta[sigma_idx]
    sigma <- exp(ls)
    mu <- as.numeric(design$X %*% beta)
    sum(.student_t_loglik(design$y, mu, sigma, nu)) +
      sum(stats::dnorm(beta, 0, 5, log = TRUE)) +
      stats::dnorm(ls, 0, 2, log = TRUE)
  }

  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)
  start <- tryCatch(
    stats::coef(stats::lm.fit(design$X, design$y)),
    error = function(e) rep.int(0, K))
  start[!is.finite(start)] <- 0

  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    init <- c(start + stats::rnorm(K, 0, 0.5), stats::rnorm(1, 0, 0.5))
    chains[[c]] <- .dram_one_chain(init, logpost, n_iter = n_iter,
                                   n_burn = n_burn,
                                   adapt_start = max(100L, n_burn %/% 2L))
  }

  draws <- .robust_t_to_draws(chains, design)
  log_lik <- .robust_t_loglik_array(chains, design, nu)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.1.0", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(koine = "0.2.0"),
    trace = list(fit = sprintf("robust-t (nu=%g) completed", nu)))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, log_likelihood = log_lik, status = "contributed")
}

# (chain, draw, var) draws with sigma back-transformed from log_sigma.
.robust_t_to_draws <- function(chains, design) {
  vnames <- c(paste0("b_", colnames(design$X)), "sigma")
  n_keep <- nrow(chains[[1L]])
  K <- design$K
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) {
    m <- chains[[c]]
    m[, K + 1L] <- exp(m[, K + 1L])
    arr[, c, ] <- m
  }
  posterior::as_draws_array(arr)
}

# Native Student-t pointwise log-likelihood as array[chain, draw, obs].
.robust_t_loglik_array <- function(chains, design, nu) {
  K <- design$K
  C <- length(chains)
  D <- nrow(chains[[1L]])
  N <- design$N
  A <- array(NA_real_, dim = c(C, D, N))
  for (c in seq_len(C)) {
    m <- chains[[c]]
    for (s in seq_len(D)) {
      mu <- as.numeric(design$X %*% m[s, seq_len(K)])
      sigma <- exp(m[s, K + 1L])
      A[c, s, ] <- .student_t_loglik(design$y, mu, sigma, nu)
    }
  }
  A
}
