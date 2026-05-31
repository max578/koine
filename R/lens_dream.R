## ----------------------------------------------------------------------------
## DREAM_ZS lens (`population / de-mcmc`) via BayesianTools
##
## Differential-Evolution Adaptive Metropolis with snooker update and past-state
## sampling (ter Braak & Vrugt 2008; Vrugt 2016). A *population* sampler: it
## evolves an ensemble of chains and proposes from their differences, a
## genuinely different geometry from the single-chain random-walk (DRAM) and the
## Gibbs lenses. It targets the same posterior (beta ~ N(0, 5); for Gaussian,
## log_sigma ~ N(0, 2)), so it is a fourth inference-family witness. It returns
## posterior draws; koine reconstructs the pointwise log-likelihood downstream
## via the family spec (F3).
## ----------------------------------------------------------------------------

#' Fit the DREAM_ZS lens
#'
#' Runs BayesianTools' DREAMzs sampler on the pass's GLM and returns an
#' [fb_fit()]. Supports the Gaussian, Binomial (logit), and Poisson (log)
#' families; for any other family/link it returns a lens-scope [abstention()].
#'
#' @param pass A `koine_pass`.
#' @param lens_spec The lens's registry record; defaults to the `dream_zs`
#'   record.
#' @param n_iter Total DREAM iterations (model evaluations across the ensemble).
#' @param burn_frac Fraction of each recorded chain discarded as burn-in.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @references ter Braak & Vrugt (2008). Differential Evolution Markov Chain
#'   with snooker updater and fewer chains. *Stat. Comput.* 18:435-446.
#' @examplesIf requireNamespace("BayesianTools", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(40))
#' d$y <- 1 + 0.8 * d$x + rnorm(40)
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit_dream(pass, n_iter = 6000L)
#' }
#' @export
fit_dream <- function(pass, lens_spec = NULL, n_iter = 15000L,
                      burn_frac = 0.5) {
  if (!requireNamespace("BayesianTools", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the BayesianTools package is not installed."))
  }
  supported <-
    (pass$family$name == "gaussian" && pass$family$link == "identity") ||
    (pass$family$name == "binomial" && pass$family$link == "logit") ||
    (pass$family$name == "poisson" && pass$family$link == "log")
  if (!supported) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("dream_zs has no setup for '%s'/'%s'.",
                       pass$family$name, pass$family$link)))
  }
  if (is.null(lens_spec)) {
    reg <- koine_stage1_registry()
    lens_spec <- reg[[which(vapply(reg, function(z) z$id, "") == "dream_zs")]]
  }
  design <- .pass_design(pass)
  has_sigma <- identical(pass$family$name, "gaussian")
  K <- design$K
  dim_theta <- K + as.integer(has_sigma)

  loglik <- function(theta) {
    beta <- theta[seq_len(K)]
    eta <- as.numeric(design$X %*% beta)
    aux <- if (has_sigma) list(sigma = exp(theta[K + 1L])) else list()
    sum(pass$family$loglik_pointwise(design$y, eta, aux))
  }
  prior_density <- function(theta) {
    lp <- sum(stats::dnorm(theta[seq_len(K)], 0, 5, log = TRUE))
    if (has_sigma) lp <- lp + stats::dnorm(theta[K + 1L], 0, 2, log = TRUE)
    lp
  }
  prior_sampler <- function(n = 1) {
    s <- cbind(matrix(stats::rnorm(n * K, 0, 5), nrow = n))
    if (has_sigma) s <- cbind(s, stats::rnorm(n, 0, 2))
    if (n == 1) as.numeric(s) else s
  }
  lower <- c(rep.int(-25, K), if (has_sigma) -8 else NULL)
  upper <- c(rep.int(25, K), if (has_sigma) 8 else NULL)

  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)
  samples <- NULL
  invisible(utils::capture.output(suppressMessages({
    prior <- BayesianTools::createPrior(
      density = prior_density, sampler = prior_sampler,
      lower = lower, upper = upper)
    setup <- BayesianTools::createBayesianSetup(
      likelihood = loglik, prior = prior, names = .dream_varnames(design,
                                                                  has_sigma))
    out <- BayesianTools::runMCMC(
      setup, sampler = "DREAMzs",
      settings = list(iterations = n_iter, message = FALSE))
    samples <- BayesianTools::getSample(out, coda = TRUE)
  })))

  draws <- .dream_to_draws(samples, design, has_sigma, burn_frac)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.0.0.9000", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(
      BayesianTools = as.character(utils::packageVersion("BayesianTools"))),
    trace = list(fit = "dream_zs completed"))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

.dream_varnames <- function(design, has_sigma) {
  v <- paste0("b_", colnames(design$X))
  if (has_sigma) c(v, "log_sigma") else v
}

# Convert DREAM's coda output (one mcmc per ensemble chain) into a draws_array,
# dropping a burn-in fraction and exponentiating log_sigma -> sigma.
.dream_to_draws <- function(samples, design, has_sigma, burn_frac) {
  chains <- if (inherits(samples, "mcmc.list")) {
    lapply(samples, function(m) as.matrix(m))
  } else {
    list(as.matrix(samples))
  }
  bnames <- paste0("b_", colnames(design$X))
  target <- if (has_sigma) c(bnames, "sigma") else bnames
  # equalise chain lengths after burn-in
  trim <- lapply(chains, function(m) {
    n <- nrow(m)
    m[seq.int(floor(n * burn_frac) + 1L, n), , drop = FALSE]
  })
  n_keep <- min(vapply(trim, nrow, integer(1)))
  arr <- array(NA_real_, dim = c(n_keep, length(trim), length(target)),
               dimnames = list(NULL, NULL, target))
  for (c in seq_along(trim)) {
    m <- trim[[c]][seq_len(n_keep), , drop = FALSE]
    src <- m[, seq_along(bnames), drop = FALSE]
    if (has_sigma) src <- cbind(src, exp(m[, design$K + 1L]))
    arr[, c, ] <- src
  }
  posterior::as_draws_array(arr)
}
