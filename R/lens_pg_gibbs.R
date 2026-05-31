## ----------------------------------------------------------------------------
## Polya-Gamma Gibbs lens (koine-authored; `gibbs / pg-exact`)
##
## An exact Gibbs sampler for Bayesian logistic regression via the Polya-Gamma
## data augmentation of Polson, Scott & Windle (2013, JASA). For y_i ~
## Bernoulli(sigmoid(x_i' beta)) with beta ~ N(0, 5^2 I):
##   omega_i | beta ~ PG(1, x_i' beta)              (drawn with BayesLogit::rpg)
##   beta | omega,y ~ N(V (X' kappa), V),  V = (X' Omega X + B0^{-1})^{-1}
##   with kappa_i = y_i - 1/2,  Omega = diag(omega).
## No tuning, no rejection --- a genuinely different inference kernel from the
## random-walk (DRAM) and Metropolis-within-Gibbs (NIMBLE) lenses. It is
## koine-authored, so it carries the `koine-implemented-validated` maturity tier
## and is discounted accordingly (SS11) until SBC-validated.
## ----------------------------------------------------------------------------

# One Polya-Gamma Gibbs chain. Returns an (n_keep x K) matrix of beta draws.
.pg_one_chain <- function(X, y, B0inv, n_iter, n_burn, beta_init) {
  N <- nrow(X)
  K <- ncol(X)
  kappa <- y - 0.5
  Xt_kappa <- crossprod(X, kappa)             # X' kappa (fixed across sweeps)
  beta <- beta_init
  keep <- matrix(NA_real_, nrow = n_iter - n_burn, ncol = K)
  for (t in seq_len(n_iter)) {
    psi <- as.numeric(X %*% beta)
    omega <- BayesLogit::rpg(N, rep.int(1, N), psi)
    V <- chol2inv(chol(crossprod(X, omega * X) + B0inv))
    m <- V %*% Xt_kappa
    beta <- as.numeric(m + crossprod(chol(V), stats::rnorm(K)))
    if (t > n_burn) keep[t - n_burn, ] <- beta
  }
  keep
}

#' Fit the Polya-Gamma Gibbs lens (logistic regression)
#'
#' Runs the exact Polya-Gamma Gibbs sampler for a Bernoulli-logit model and
#' returns an [fb_fit()]. Resolves only for the binomial (logit) family; for any
#' other family it returns a lens-scope [abstention()] rather than a wrong fit.
#'
#' @param pass A `koine_pass` with a binomial (logit) family.
#' @param lens_spec The lens's registry record; defaults to the `pg_gibbs`
#'   record.
#' @param n_iter Iterations per chain.
#' @param n_burn Burn-in iterations discarded per chain.
#' @param n_chains Number of chains.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @references Polson, Scott & Windle (2013). Bayesian inference for logistic
#'   models using Polya-Gamma latent variables. *JASA* 108(504):1339-1349.
#' @examplesIf requireNamespace("BayesLogit", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(80))
#' d$y <- rbinom(80, 1, plogis(0.5 * d$x))
#' pass <- koine_pass(y ~ x, d, binomial())
#' fit_pg_gibbs(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_pg_gibbs <- function(pass, lens_spec = NULL, n_iter = 3000L,
                         n_burn = 1500L, n_chains = 4L) {
  if (!requireNamespace("BayesLogit", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the BayesLogit package is not installed."))
  }
  if (!(pass$family$name == "binomial" && pass$family$link == "logit")) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("pg_gibbs serves binomial/logit only, not '%s'/'%s'.",
                       pass$family$name, pass$family$link)))
  }
  if (is.null(lens_spec)) {
    reg <- koine_stage1_registry()
    lens_spec <- reg[[which(vapply(reg, function(z) z$id, "") == "pg_gibbs")]]
  }
  design <- .pass_design(pass)
  if (!all(design$y %in% c(0, 1))) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = "pg_gibbs supports Bernoulli (0/1) outcomes only."))
  }
  B0inv <- diag(1 / 25, design$K)             # beta ~ N(0, 5^2): B0 = 25 I
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  set.seed(dseed)
  warm <- tryCatch(
    suppressWarnings(stats::coef(stats::glm.fit(
      design$X, design$y, family = stats::binomial()))),
    error = function(e) rep.int(0, design$K))
  warm[!is.finite(warm)] <- 0

  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    init <- warm + stats::rnorm(design$K, 0, 0.5)
    chains[[c]] <- .pg_one_chain(design$X, design$y, B0inv,
                                 n_iter = n_iter, n_burn = n_burn,
                                 beta_init = init)
  }
  draws <- .beta_chains_to_draws(chains, design)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.0.0.9000", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(
      BayesLogit = as.character(utils::packageVersion("BayesLogit"))),
    trace = list(fit = "polya-gamma gibbs completed"))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

# Coefficient-only chains (no aux param) -> posterior draws_array with b_<col>
# variable names.
.beta_chains_to_draws <- function(chains, design) {
  vnames <- paste0("b_", colnames(design$X))
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) arr[, c, ] <- chains[[c]]
  posterior::as_draws_array(arr)
}
