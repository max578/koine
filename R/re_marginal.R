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

# Weakly-informative prior SDs for the RE lenses, per family. The shared
# N(0, 5) coefficient prior is fine for the Gaussian (identity) and binomial
# (logit, bounded outcome) families, but for the Poisson (log link, unbounded
# count) it makes the prior-predictive pathological (simulated counts ~1e22), so
# SBC cannot validate the lens. A tighter family-appropriate prior gives a
# sensible count prior-predictive. Both RE lenses share these SDs per family, so
# their agreement remains method-independent (the priors match each other).
.re_prior <- function(family_name) {
  if (identical(family_name, "poisson")) {
    list(beta_sd = 1, log_sigma_sd = 2, log_tau_sd = 0.5)
  } else {
    list(beta_sd = 5, log_sigma_sd = 2, log_tau_sd = 2)
  }
}

# Whether the Gaussian-LMM analytic shortcut applies (identity-link Gaussian:
# the random intercept enters the mean linearly, so the marginal is exactly
# normal). Any other family/link routes through Gauss-Hermite quadrature.
.re_is_gaussian_lmm <- function(family_spec) {
  identical(family_spec$name, "gaussian") &&
    identical(family_spec$link, "identity")
}

# Adaptive Gauss-Hermite quadrature of a random-effect integral, per "unit"
# (a unit is a group of observations sharing one random effect b_u; for the
# per-observation predictive each observation is its own unit). For each unit it
# Newton-iterates to the mode b-hat_u of f_u(b) = sum_{i in u} loglik(y_i,
# eta_i + b) + log N(b | 0, tau^2) using the family's canonical-link score /
# information, scales standard Gauss-Hermite nodes by the local curvature, and
# returns the per-unit log marginal log integral exp(f_u(b)) db. FIXED Gauss-
# Hermite is non-convergent for sharply-peaked count integrands (the Poisson
# group marginal); centring the nodes on the integrand makes the quadrature
# accurate with few nodes. The quadrature still converges to the exact integral
# --- it is a node placement, not a Laplace approximation to the posterior.
.aghq_marginal <- function(family_spec, y, eta_fixed, unit, tau, aux,
                           gh = NULL, return_post = FALSE) {
  if (is.null(gh)) {
    gh <- .gh_rule(20L)
  }
  uf <- as.integer(factor(unit))
  n_u <- max(uf)
  nq <- length(gh$nodes)
  inv_t2 <- 1 / tau^2
  b <- numeric(n_u)
  # Newton-iterate to the per-unit mode b-hat. Run to a tight tolerance: a
  # stationary mode makes the node placement b-hat(theta) / sd-hat(theta) locally
  # insensitive to theta, so the score-identity gradient below (which holds the
  # nodes fixed) matches the gradient of the computed quadrature to quadrature
  # order --- the condition the finite-difference gradient gate verifies.
  for (it in seq_len(50L)) {
    eta <- eta_fixed + b[uf]
    d1 <- family_spec$dlp_deta(y, eta, aux)
    d2 <- family_spec$d2lp_deta2(y, eta, aux)
    # extreme theta along a Hamiltonian trajectory can overflow the GLM mean
    # (e.g. exp(eta) for the log-link Poisson). Bail rather than crash: the
    # downstream value is then non-finite and the lens marks a divergence.
    if (!all(is.finite(d1)) || !all(is.finite(d2))) break
    g1 <- as.numeric(rowsum(d1, uf)) - b * inv_t2
    g2 <- as.numeric(rowsum(d2, uf)) - inv_t2
    step <- g1 / g2
    if (!all(is.finite(step))) break
    b <- b - step
    if (max(abs(step)) < 1e-10) break
  }
  info <- -(as.numeric(rowsum(family_spec$d2lp_deta2(y, eta_fixed + b[uf], aux),
                              uf)) - inv_t2)
  sd_hat <- sqrt(1 / info)
  log_w <- log(gh$weights)
  # log integrand at every (unit, node): the n_u x nq matrix of fu + log weight.
  b_nodes <- matrix(NA_real_, nrow = n_u, ncol = nq)
  lcols <- matrix(NA_real_, nrow = n_u, ncol = nq)
  for (k in seq_len(nq)) {
    bk <- b + sqrt(2) * sd_hat * gh$nodes[k]
    b_nodes[, k] <- bk
    fu <- as.numeric(rowsum(family_spec$loglik_pointwise(y, eta_fixed + bk[uf],
                                                         aux), uf)) +
      stats::dnorm(bk, 0, tau, log = TRUE)
    lcols[, k] <- fu + log_w[k] + gh$nodes[k]^2
  }
  lse <- .row_log_sum_exp(lcols)
  value <- log(sqrt(2) * sd_hat) + lse
  if (!return_post) {
    return(value)
  }
  # posterior weights over the quadrature nodes (the discrete tilted density
  # q_u, normalised per unit): w_post[u, k] proportional to exp(lcols[u, k]).
  w_post <- exp(lcols - lse)
  list(value = value, uf = uf, b_nodes = b_nodes, w_post = w_post)
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
  # each observation is its own unit (a fresh random effect, marginalised)
  .aghq_marginal(family_spec, y, eta_fixed, seq_along(y), tau, aux, gh)
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
  # adaptive Gauss-Hermite per group (fixed GH is non-convergent for the
  # sharply-peaked Poisson group integrand)
  sum(.aghq_marginal(family_spec, y, eta_fixed, group, tau, aux, gh))
}

# Exact Gaussian-LMM per-group marginal: y_g ~ N(eta_fixed_g, sigma^2 I + tau^2
# 11'). Uses the matrix-determinant lemma and Sherman-Morrison so no per-group
# matrix is formed: for a group of size m with residual r = y_g - eta_fixed_g,
#   log|Sigma| = m log(sigma^2) + log(1 + m tau^2 / sigma^2)
#   r'Sigma^-1 r = (sum r^2 - tau^2/(sigma^2 + m tau^2) (sum r)^2) / sigma^2.
# Per-group Gaussian-LMM log marginal (a vector over groups). The scalar version
# below sums it; the t random-effect mixture (below) needs the per-group values
# because each group's random effect mixes over its own scale.
.re_group_gaussian_vec <- function(y, eta_fixed, group, tau, sigma) {
  s2 <- sigma^2
  t2 <- tau^2
  r <- y - eta_fixed
  m <- as.numeric(table(group))
  sum_r <- as.numeric(rowsum(r, group))
  sum_r2 <- as.numeric(rowsum(r^2, group))
  log_det <- m * log(s2) + log1p(m * t2 / s2)
  quad <- (sum_r2 - (t2 / (s2 + m * t2)) * sum_r^2) / s2
  -0.5 * (m * log(2 * pi) + log_det + quad)
}

.re_group_marginal_gaussian <- function(y, eta_fixed, group, tau, sigma) {
  sum(.re_group_gaussian_vec(y, eta_fixed, group, tau, sigma))
}

## ----------------------------------------------------------------------------
## Gradient of the per-group marginal log-likelihood (Stage-2; HMC kernel)
##
## The Hamiltonian marginal lens needs the gradient of the marginal log-
## likelihood with respect to the fixed-part linear predictor (so the chain rule
## delivers d/d beta = X' s) and with respect to log tau (and log sigma for the
## Gaussian LMM). Two routes, matching the two value routes:
##   * Gaussian LMM --- the marginal is exactly normal per group, so the gradient
##     is closed form (the GLS residual gives s; the variance-component
##     derivatives come from the matrix-determinant / Sherman-Morrison identities
##     used for the value).
##   * GLMM --- the score (Fisher) identity: d/d theta log int exp(f(b)) db =
##     E_{q(b)}[d/d theta f(b)], with q(b) propto exp(f(b)) the tilted density.
##     The adaptive Gauss-Hermite nodes and their normalised posterior weights
##     (`return_post`) approximate that expectation directly, so one quadrature
##     pass returns both the value and the gradient. Verified against finite
##     differences of the value (the `test-re_marginal` gradient gate).
## Returns: `value` (scalar), `s_obs` (length-N expected score per observation),
## `d_logtau` (d value / d log tau, holding eta_fixed fixed), and `d_logsigma`
## (Gaussian only; `NA` otherwise).
## ----------------------------------------------------------------------------
.re_group_marginal_grad <- function(family_spec, y, eta_fixed, group, tau, aux,
                                    gh = NULL) {
  if (.re_is_gaussian_lmm(family_spec)) {
    return(.re_group_gaussian_grad(y, eta_fixed, group, tau, aux$sigma))
  }
  post <- .aghq_marginal(family_spec, y, eta_fixed, group, tau, aux, gh,
                         return_post = TRUE)
  uf <- post$uf
  n_u <- nrow(post$w_post)
  nq <- ncol(post$w_post)
  s_obs <- numeric(length(y))
  eb2 <- numeric(n_u)
  for (k in seq_len(nq)) {
    bk <- post$b_nodes[, k]
    w_obs <- post$w_post[uf, k]
    s_obs <- s_obs +
      w_obs * family_spec$dlp_deta(y, eta_fixed + bk[uf], aux)
    eb2 <- eb2 + post$w_post[, k] * bk^2
  }
  list(value = sum(post$value), s_obs = s_obs,
       d_logtau = sum(eb2 / tau^2 - 1), d_logsigma = NA_real_)
}

# Closed-form gradient of the exact Gaussian-LMM per-group marginal. `s_obs` is
# the GLS residual Sigma_g^-1 (y_g - eta_g) per observation (so d/d beta = X' s);
# the variance-component derivatives are d value / d log sigma and d value / d
# log tau via the compound-symmetry determinant / quadratic-form derivatives.
.re_group_gaussian_grad <- function(y, eta_fixed, group, tau, sigma) {
  s2 <- sigma^2
  t2 <- tau^2
  r <- y - eta_fixed
  gf <- as.integer(factor(group))
  m <- as.numeric(table(gf))                       # group sizes
  sum_r <- as.numeric(rowsum(r, gf))
  d_grp <- s2 + m * t2                             # D_g = sigma^2 + m_g tau^2
  # GLS residual per observation: (1/s2)(r_i - (t2 / D_g) sum_r_g)
  s_obs <- (r - (t2 / d_grp)[gf] * sum_r[gf]) / s2
  sum_r2 <- as.numeric(rowsum(r^2, gf))
  b_quad <- sum_r^2                                # B_g = (sum r_g)^2
  # d log p_g / d t2 = -0.5 (m/D) + 0.5 B/D^2 ; d/d log tau = sum_g (.) * 2 t2
  dlp_dt2 <- -0.5 * m / d_grp + 0.5 * b_quad / d_grp^2
  d_logtau <- sum(dlp_dt2) * 2 * t2
  # d log p_g / d s2, then d/d log sigma = sum_g (.) * 2 s2
  dlogdet_ds2 <- m / s2 - m * t2 / (s2 * d_grp)
  dquad_ds2 <- -sum_r2 / s2^2 +
    t2 * b_quad * (d_grp + s2) / (s2 * d_grp)^2
  dlp_ds2 <- -0.5 * (dlogdet_ds2 + dquad_ds2)
  d_logsigma <- sum(dlp_ds2) * 2 * s2
  list(value = .re_group_marginal_gaussian(y, eta_fixed, group, tau, sigma),
       s_obs = s_obs, d_logtau = d_logtau, d_logsigma = d_logsigma)
}

# Generalized Gauss-Laguerre rule (weight x^alpha exp(-x) on [0, Inf)) of order
# `n`, by the Golub-Welsch eigendecomposition of the Jacobi matrix. Integrates
# integral g(x) x^alpha exp(-x) dx ~= sum_k weights_k g(nodes_k), with
# sum_k weights_k = gamma(alpha + 1).
.gl_rule <- function(n = 32L, alpha = 0) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(nodes = alpha + 1, weights = gamma(alpha + 1)))
  }
  k <- seq_len(n)
  diagv <- 2 * (k - 1L) + alpha + 1
  i <- seq_len(n - 1L)
  off <- sqrt(i * (i + alpha))
  jac <- matrix(0, n, n)
  diag(jac) <- diagv
  jac[cbind(i, i + 1L)] <- off
  jac[cbind(i + 1L, i)] <- off
  e <- eigen(jac, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord],
       weights = (gamma(alpha + 1) * e$vectors[1L, ]^2)[ord])
}

# Quadrature for the Student-t scale mixture: a t_nu random effect is N(0,
# tau^2 / lambda) with lambda ~ Gamma(nu/2, nu/2). Returns the lambda nodes and
# normalised weights so that E_lambda[h(lambda)] ~= sum_k weight_k h(lambda_k).
.t_mix_rule <- function(nu, n = 32L) {
  a <- nu / 2
  gl <- .gl_rule(n, alpha = a - 1)
  list(lambda = gl$nodes / a, weight = gl$weights / gamma(a))
}

# Per-observation integrated predictive for a Gaussian LMM with t_nu random
# intercepts: log E_lambda[ N(y_i | eta_fixed_i, sigma^2 + tau^2 / lambda) ],
# the honest predictive for a new observation whose heavy-tailed random effect
# is marginalised out.
.re_integrated_loglik_t <- function(y, eta_fixed, tau, sigma, nu, mix = NULL) {
  if (is.null(mix)) {
    mix <- .t_mix_rule(nu)
  }
  cols <- vapply(seq_along(mix$lambda), function(k) {
    sd_k <- sqrt(sigma^2 + tau^2 / mix$lambda[k])
    stats::dnorm(y, eta_fixed, sd_k, log = TRUE) + log(mix$weight[k])
  }, numeric(length(y)))
  .row_log_sum_exp(matrix(cols, nrow = length(y)))
}

# Per-group marginal log-likelihood for a Gaussian LMM with t_nu random
# intercepts. Each group's random effect mixes over its OWN scale, so the
# integral is per group: sum_g log E_lambda[ N(y_g | eta_g, sigma^2 I + (tau^2 /
# lambda) 11') ], with the inner Gaussian marginal evaluated analytically.
.re_group_marginal_t <- function(y, eta_fixed, group, tau, sigma, nu,
                                 mix = NULL) {
  if (is.null(mix)) {
    mix <- .t_mix_rule(nu)
  }
  ng <- length(unique(group))
  per_group <- vapply(seq_along(mix$lambda), function(k) {
    tau_k <- tau / sqrt(mix$lambda[k])
    .re_group_gaussian_vec(y, eta_fixed, group, tau_k, sigma) +
      log(mix$weight[k])
  }, numeric(ng))
  sum(.row_log_sum_exp(matrix(per_group, nrow = ng)))
}
