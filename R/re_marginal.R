## ----------------------------------------------------------------------------
## Random-effects marginalisation --- the shared integration core (Stage-2)
##
## A single random intercept enters the linear predictor as eta_i = x_i' beta +
## b_{g(i)}, with b_g ~ N(0, tau^2). koine never samples b when it forms the
## fourth opinion: it integrates b out. Two distinct integrals are needed, both
## served from this one place so the lens and the predictive can never drift
## apart:
##
##   * the per-GROUP marginal likelihood (lens fitting) integrates b_g once over
##     the *product* of the group's observations --- exact multivariate-normal
##     for the Gaussian LMM, Gauss-Hermite quadrature for a GLMM;
##   * the per-OBSERVATION integrated predictive (LOO stacking) integrates b for
##     a *single* held-out observation --- analytic N(y_i | x_i'beta, sigma^2 +
##     tau^2) for the Gaussian LMM, the same quadrature otherwise.
##
## Standardising every RE lens onto the per-observation integrated predictive is
## what makes the stacking weights valid: all components then predict the same
## held-out quantity (a new observation, its random effect marginalised under
## the sampled variance components), not one conditioned on a random effect the
## held-out point helped estimate (which would inflate the LOO and push the PSIS
## Pareto-k past the reliability gate).
## ----------------------------------------------------------------------------

# Gauss-Hermite rule (physicists', weight exp(-x^2)) of order `n`, by the
# Golub-Welsch eigendecomposition of the symmetric tridiagonal Jacobi matrix.
# Returns nodes and weights for the approximation integral f(x) exp(-x^2) dx
# ~= sum_k weights_k f(nodes_k), with sum_k weights_k = sqrt(pi).
.gh_rule <- function(n = 20L) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(nodes = 0, weights = sqrt(pi)))
  }
  i <- seq_len(n - 1L)
  off <- sqrt(i / 2)
  jac <- matrix(0, n, n)
  jac[cbind(i, i + 1L)] <- off
  jac[cbind(i + 1L, i)] <- off
  e <- eigen(jac, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = (sqrt(pi) * e$vectors[1L, ]^2)[ord])
}

# Per-row log-sum-exp of a matrix (rows = observations or groups, cols = nodes).
.row_log_sum_exp <- function(mat) {
  m <- apply(mat, 1L, max)
  m + log(rowSums(exp(mat - m)))
}

# Whether the Gaussian-LMM analytic shortcut applies (identity-link Gaussian:
# the random intercept enters the mean linearly, so the marginal is exactly
# normal). Any other family/link routes through Gauss-Hermite quadrature.
.re_is_gaussian_lmm <- function(family_spec) {
  identical(family_spec$name, "gaussian") &&
    identical(family_spec$link, "identity")
}

#' Per-observation integrated random-effect predictive (Stage-2)
#'
#' Computes `log integral p(y_i | linkinv(eta_fixed_i + b)) N(b | 0, tau^2) db`
#' for every observation --- the honest predictive density for a *new*
#' observation with its random effect marginalised out. Analytic for the
#' identity-link Gaussian LMM; Gauss-Hermite quadrature otherwise. This is the
#' quantity every RE lens is stacked on, so the LOO weights are comparable
#' across lenses (the same held-out target).
#'
#' @param family_spec A `koine_family_spec`.
#' @param y Response vector (length `N`).
#' @param eta_fixed Fixed-part linear predictor `x_i' beta` (length `N`).
#' @param tau Random-intercept standard deviation (scalar, `>= 0`).
#' @param aux Named list of auxiliary parameters (e.g. `list(sigma = .)`).
#' @param gh A Gauss-Hermite rule from `.gh_rule()`; computed if `NULL`.
#' @return A length-`N` vector of integrated pointwise log-likelihoods.
#' @keywords internal
.re_integrated_loglik <- function(family_spec, y, eta_fixed, tau, aux,
                                  gh = NULL) {
  if (.re_is_gaussian_lmm(family_spec)) {
    sd_marg <- sqrt(aux$sigma^2 + tau^2)
    return(stats::dnorm(y, mean = eta_fixed, sd = sd_marg, log = TRUE))
  }
  if (is.null(gh)) {
    gh <- .gh_rule()
  }
  shift <- sqrt(2) * tau
  log_w <- log(gh$weights / sqrt(pi))
  cols <- vapply(seq_along(gh$nodes), function(k) {
    family_spec$loglik_pointwise(y, eta_fixed + shift * gh$nodes[k], aux) +
      log_w[k]
  }, numeric(length(y)))
  .row_log_sum_exp(matrix(cols, nrow = length(y)))
}

#' Per-group random-effect marginal log-likelihood (Stage-2)
#'
#' Computes `sum_g log int_b [prod_{i in g} p(y_i|mu_i(b))] N(b|0,tau^2) db`
#' (`mu_i(b) = linkinv(eta_fixed_i + b)`) --- the marginal likelihood the RE
#' lenses target when they
#' sample `(beta, sigma, tau)` without ever sampling the random effects. Exact
#' for the identity-link Gaussian LMM (a multivariate normal per group, via the
#' determinant / Sherman-Morrison identities); Gauss-Hermite quadrature
#' otherwise.
#'
#' @param family_spec A `koine_family_spec`.
#' @param y Response vector (length `N`).
#' @param eta_fixed Fixed-part linear predictor (length `N`).
#' @param group Integer group index per observation (length `N`).
#' @param tau Random-intercept standard deviation (scalar, `>= 0`).
#' @param aux Named list of auxiliary parameters (e.g. `list(sigma = .)`).
#' @param gh A Gauss-Hermite rule from `.gh_rule()`; computed if `NULL`.
#' @return A scalar log marginal likelihood.
#' @keywords internal
.re_group_marginal_loglik <- function(family_spec, y, eta_fixed, group, tau,
                                      aux, gh = NULL) {
  if (.re_is_gaussian_lmm(family_spec)) {
    return(.re_group_marginal_gaussian(y, eta_fixed, group, tau, aux$sigma))
  }
  if (is.null(gh)) {
    gh <- .gh_rule()
  }
  shift <- sqrt(2) * tau
  log_w <- log(gh$weights / sqrt(pi))
  # For each node, the per-group sum of pointwise log-likelihoods (rowsum groups
  # the per-observation contributions); columns are then the integrand summands.
  n_groups <- length(unique(group))
  group_sums <- vapply(seq_along(gh$nodes), function(k) {
    pw <- family_spec$loglik_pointwise(y, eta_fixed + shift * gh$nodes[k], aux)
    as.numeric(rowsum(pw, group)) + log_w[k]
  }, numeric(n_groups))
  sum(.row_log_sum_exp(matrix(group_sums, nrow = n_groups)))
}

# Exact Gaussian-LMM per-group marginal: y_g ~ N(eta_fixed_g, sigma^2 I + tau^2
# 11'). Uses the matrix-determinant lemma and Sherman-Morrison so no per-group
# matrix is formed: for a group of size m with residual r = y_g - eta_fixed_g,
#   log|Sigma| = m log(sigma^2) + log(1 + m tau^2 / sigma^2)
#   r'Sigma^-1 r = (sum r^2 - tau^2/(sigma^2 + m tau^2) (sum r)^2) / sigma^2.
.re_group_marginal_gaussian <- function(y, eta_fixed, group, tau, sigma) {
  s2 <- sigma^2
  t2 <- tau^2
  r <- y - eta_fixed
  m <- as.numeric(table(group))
  sum_r <- as.numeric(rowsum(r, group))
  sum_r2 <- as.numeric(rowsum(r^2, group))
  log_det <- m * log(s2) + log1p(m * t2 / s2)
  quad <- (sum_r2 - (t2 / (s2 + m * t2)) * sum_r^2) / s2
  -0.5 * sum(m * log(2 * pi) + log_det + quad)
}
