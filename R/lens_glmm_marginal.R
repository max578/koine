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
## The Gaussian and binomial marginals ride the gradient-informed Hamiltonian
## (NUTS) kernel. The Poisson marginal's variance-component funnel defeats a
## static-metric Hamiltonian kernel, so Poisson rides a funnel-safe BLOCKED
## double-slice kernel instead (slice log tau; slice the coefficients along their
## tau-resolved Laplace principal axes); its unmixable large-tau tail is caught by
## a post-fit convergence gate that makes the lens abstain (see the blocked-kernel
## and `fit_glmm_marginal` notes below). Either way the lens *quadrature-integrates*
## the random-effect integral --- deliberately orthogonal to the sealed engines:
## NIMBLE *samples* the random effects (Gibbs), INLA *Laplace-approximates* the
## integral, and this lens *quadrature-integrates* it --- three distinct treatments
## of the one integral, so their agreement (or disagreement) carries a genuine
## fourth-witness signal, strongest for GLMMs where the Laplace approximation is
## only approximate. It reports `(beta, sigma, tau)` draws and leaves the pointwise
## log-likelihood for `assemble_loglik` to compute as the per-observation
## integrated predictive (F3), so it is stacked on the identical held-out target as
## the NIMBLE GLMM lens.
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

# Value-ONLY marginal log-posterior (the scalar the blocked slice kernel needs).
# It computes the marginal log-likelihood directly (`.re_group_marginal_loglik`)
# and skips the score-identity gradient quadrature entirely --- roughly half the
# cost of the value-and-gradient path, and the slice sampler never uses the
# gradient. The arithmetic is identical to the `value` branch of
# `.glmm_marginal_logpost_grad` (same reparameterisation, same priors), so the two
# agree to numerical precision (asserted in `test-glmm_marginal`).
.glmm_marginal_logpost_value <- function(theta, design, redesign, pass, lay, gh,
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
  ll <- .re_group_marginal_loglik(pass$family, design$y, eta_fixed,
                                  redesign$group, tau, aux, gh)
  value <- ll + sum(stats::dnorm(beta, 0, prior$beta_sd, log = TRUE)) + lp_extra
  if (!is.finite(value)) {
    return(-Inf)
  }
  value
}

## ----------------------------------------------------------------------------
## Blocked double-slice kernel for the Poisson marginal (funnel-safe)
##
## The Poisson marginal posterior over (beta, log tau) is a variance-component
## FUNNEL: the conditioning of the geometry varies by position (the AGHQ-probed
## condition number runs ~7 -> 124 across the log-tau range, with the eigenframe
## stable --- a position-dependent SCALING, not a rotation). A single static-metric
## NUTS cannot serve the neck and the body at once: it collapses the global step
## size, giving ~50% divergences, ess_tau ~ 3, and order-minutes (capped) fits on
## the minority of datasets with a large random-intercept SD. Marginalising the
## random effects rules out the usual non-centred remedy (there is no b to
## reparameterise).
##
## The fix BLOCKS the funnel coordinate off from a gradient kernel and slices both
## blocks:
##   * log tau is updated by a 1-D slice sampler (`.slice_1d`) --- no step size to
##     mistune, so it brackets the conditional and moves at the correct scale at
##     every point of the funnel;
##   * the coefficients are slice-sampled at FIXED tau along their tau-resolved
##     principal axes (the scaled eigenvectors of the conditional's Laplace
##     covariance), where the soft funnel direction vanishes and the conditional is
##     well-conditioned. A gradient (NUTS) coefficient update funnels into depth-
##     saturated trees here; a Laplace independence proposal freezes (acceptance
##     collapses); a slice along the tau-resolved axes cannot stall and needs no
##     tuning.
## The composition leaves the marginal posterior invariant (each block is a valid
## tau-/beta-conditional slice transition). The level reparameterisation
## (beta0 = alpha0 - tau^2/2) composes cleanly: slicing log tau with the level
## coordinate alpha0 held fixed is exactly the intercept/scale decorrelation a
## Gibbs block wants. It keeps the lens's core independence (the marginalisation +
## quadrature treatment of the RE integral) intact; only the outer kernel changes,
## and only for Poisson. The large-tau tail that no marginal sampler mixes is
## caught by `fit_glmm_marginal`'s post-fit convergence gate (the lens abstains).
## ----------------------------------------------------------------------------

# The coefficient-block target at a FIXED log tau: theta -> list(value, grad over
# the coefficient sub-vector), assembled by pinning the tau coordinate of the
# full marginal log-posterior. The gradient path anchors the tau-knot Laplace
# axes (the conditional mode + curvature); the slice itself rides the value-only
# conditional.
.glmm_marginal_alpha_lpg <- function(lpg_full, lay, lt) {
  function(alpha) {
    theta <- numeric(lay$dim)
    theta[lay$coef_idx] <- alpha
    theta[lay$tau_idx] <- lt
    vg <- lpg_full(theta)
    list(value = vg$value, grad = vg$grad[lay$coef_idx])
  }
}

# Value-only coefficient conditional at a FIXED log tau (the scalar the
# coefficient slice rides), assembled from the value-only marginal log-posterior.
.glmm_alpha_value_fn <- function(lpg_val, lay, lt) {
  function(alpha) {
    theta <- numeric(lay$dim)
    theta[lay$coef_idx] <- alpha
    theta[lay$tau_idx] <- lt
    lpg_val(theta)
  }
}

# Symmetric solve with a regularised fallback: returns `A^-1` (or the identity if
# `A` is singular / non-finite), so a degenerate conditional curvature never
# crashes the metric build.
.sym_solve <- function(a_mat, k) {
  out <- tryCatch(solve(a_mat), error = function(e) NULL)
  if (is.null(out) || !all(is.finite(out))) {
    return(diag(k))
  }
  (out + t(out)) / 2
}

# Finite-difference curvature (precision) of the coefficient conditional: the
# negative Hessian of the conditional log-posterior value, by central differences
# of the analytic coefficient gradient. Symmetrised; this is the local precision
# whose inverse is the coefficient-block mass matrix.
.glmm_alpha_hess <- function(la, a, h = 1e-4) {
  k <- length(a)
  hmat <- matrix(0, k, k)
  for (j in seq_len(k)) {
    ap <- a; ap[j] <- ap[j] + h
    am <- a; am[j] <- am[j] - h
    hmat[, j] <- (la(ap)$grad - la(am)$grad) / (2 * h)
  }
  -(hmat + t(hmat)) / 2
}

# Newton mode of the coefficient conditional, warm-started from `a`. A handful of
# damped Newton steps (backtracking if a step leaves the finite-value region)
# locate the mode; the conditional is smooth and unimodal, so this converges in a
# few steps. Used only to anchor the tau-adaptive metric at a tau-only reference.
.glmm_alpha_mode <- function(la, a, newton = 6L) {
  for (it in seq_len(newton)) {
    g <- la(a)$grad
    if (!all(is.finite(g))) break
    step <- tryCatch(solve(.glmm_alpha_hess(la, a), g),
                     error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
    a_new <- a + step
    if (!is.finite(la(a_new)$value)) {
      a_new <- a + 0.25 * step                   # backtrack toward the interior
      if (!is.finite(la(a_new)$value)) break
    }
    a <- a_new
    if (max(abs(step)) < 1e-4) break
  }
  a
}

## ----------------------------------------------------------------------------
## tau-adaptive coefficient slice axes (the funnel-safe, cost-flat blocked kernel)
##
## Blocking the variance component (slice log tau) from the coefficients (slice
## them at fixed tau) clears the tau ~ 1.5 neck where a single static-metric NUTS
## funnels (~50% divergences, ess_tau ~ 3) and a Laplace independence proposal
## freezes (acceptance collapses). koine slice-samples log tau, then slice-samples
## the coefficients along the tau-resolved principal axes of their conditional's
## Laplace covariance. The Laplace approximation supplies only the AXES (its scaled
## eigenvectors), which decorrelate and rescale the block so each 1-D slice has
## unit width; the slice itself --- no proposal, no step size --- brackets the true
## conditional and cannot stall even where the Gaussian approximation is poor.
## Correctness never depends on the approximation being good, only on the axes
## spanning the block; the target is the EXACT-quadrature marginal, so the lens's
## core independence (exact AGHQ treatment of the random-effect integral, distinct
## from INLA's Laplace approximation) is preserved --- the approximation only
## orients the slice.
##
## At a large random-intercept SD (tau beyond ~2, the ~8% upper tail of the prior)
## the marginalised funnel's soft direction is an (intercept, log tau) diagonal
## that no marginal sampler traverses cheaply --- with the random effects integrated
## out, the usual non-centred remedy is unavailable. The lens does not pretend
## otherwise: the fit then under-mixes the variance component, and a post-fit
## convergence gate (`fit_glmm_marginal`) makes the lens ABSTAIN rather than
## contribute a badly-mixed opinion (a first-class, honest outcome, SS3). koine
## contributes a calibrated Poisson opinion where the variance component is
## identifiable, and abstains where it is not.
##
## The axes are built once (during warm-up) at a fixed grid of log-tau knots
## spanning the prior-plausible range, each anchored at that knot's conditional
## mode (knots walked in order, each warm-started from its neighbour); sampling
## looks up the nearest knot in O(1).
## ----------------------------------------------------------------------------
.glmm_build_tau_knots <- function(lpg_full, lay, alpha0,
                                  lt_range = c(log(0.15), log(7)),
                                  n_knots = 9L, inflate = 2) {
  ka <- length(lay$coef_idx)
  knot_lt <- seq(lt_range[1L], lt_range[2L], length.out = n_knots)
  axes <- vector("list", n_knots)
  a <- alpha0
  for (i in seq_len(n_knots)) {
    la <- .glmm_marginal_alpha_lpg(lpg_full, lay, knot_lt[i])
    a <- .glmm_alpha_mode(la, a)                 # warm-start from the previous knot
    prec <- .glmm_alpha_hess(la, a)
    cov <- .nuts_regularise(.sym_solve(prec, ka), ka)
    e <- eigen(cov, symmetric = TRUE)
    # scaled principal axes (columns): a unit step along axis j moves `inflate`
    # local standard deviations in that conditional direction. The mild over-width
    # (`inflate`) keeps the slice's stepping-out cheap when the large-tau
    # conditional is heavier-tailed than its Laplace approximation (over-wide costs
    # logarithmic shrink steps, under-wide costs linear step-outs). Eigenvalues are
    # floored AND capped so a degenerate mode estimate never yields a zero / runaway
    # slice width.
    sds <- inflate * sqrt(pmin(pmax(e$values, 1e-6), 1e4))
    axes[[i]] <- e$vectors %*% diag(sds, ka)
  }
  list(lt = knot_lt, axes = axes)
}

# Nearest-knot coefficient slice axes for the current log tau (a step function of
# tau, hence a valid tau-only orientation).
.glmm_knot_axes <- function(knots, lt) {
  knots$axes[[which.min(abs(knots$lt - lt))]]
}

# Slice-sample the coefficient block along its tau-resolved principal axes (random
# axis order), each a 1-D slice through the current point with unit width (the
# axis scale is carried in `axes`). `la_val` is the value-only coefficient
# conditional. Leaves the conditional invariant and cannot get stuck.
.glmm_slice_alpha <- function(la_val, alpha, axes) {
  ka <- length(alpha)
  for (j in sample.int(ka)) {
    u <- axes[, j]
    along <- function(s) la_val(alpha + s * u)
    s <- .slice_1d(along, 0, w = 1)
    alpha <- alpha + s * u
  }
  alpha
}

# One blocked chain: slice log tau, then slice the coefficients along their
# tau-resolved Laplace principal axes. The axes (`knots`) are built from this
# chain's warm start (`lpg_full` supplies the conditional mode + curvature); the
# log-tau slice width is set from the warm-up spread. Returns the (kept x dim) draw
# matrix in the same layout the joint kernel returns (coefficients on the sampled
# level scale, log tau on the log scale).
.glmm_marginal_blocked_chain <- function(init, lpg_full, lpg_val, lay, n_iter,
                                         n_warmup) {
  alpha <- init[lay$coef_idx]
  lt <- init[lay$tau_idx]
  # the log-tau slice rides the value-only marginal (the gradient is never used).
  tau_slice <- function(a) function(z) {
    theta <- numeric(lay$dim)
    theta[lay$coef_idx] <- a
    theta[lay$tau_idx] <- z
    lpg_val(theta)
  }
  knots <- .glmm_build_tau_knots(lpg_full, lay, alpha)

  # one sweep: slice log tau | alpha, then slice the coefficients | tau along their
  # tau-resolved principal axes. Both are slice updates on the value-only marginal,
  # so the sweep needs no step-size or acceptance tuning.
  sweep <- function(alpha, lt, w) {
    lt <- .slice_1d(tau_slice(alpha), lt, w = w, lower = -Inf, upper = Inf)
    la_val <- .glmm_alpha_value_fn(lpg_val, lay, lt)
    alpha <- .glmm_slice_alpha(la_val, alpha, .glmm_knot_axes(knots, lt))
    list(alpha = alpha, lt = lt)
  }

  n_keep <- n_iter - n_warmup
  keep <- matrix(NA_real_, nrow = n_keep, ncol = lay$dim)
  w <- 1
  seen_lt <- numeric(n_warmup)
  for (m in seq_len(n_warmup)) {
    s <- sweep(alpha, lt, w)
    alpha <- s$alpha
    lt <- s$lt
    seen_lt[m] <- lt
    if (m == n_warmup %/% 2L) {                  # adapt the log-tau slice width once
      sw <- stats::sd(seen_lt[seq_len(m)][is.finite(seen_lt[seq_len(m)])])
      if (is.finite(sw) && sw > 0) w <- sw
    }
  }
  for (i in seq_len(n_keep)) {
    s <- sweep(alpha, lt, w)
    alpha <- s$alpha
    lt <- s$lt
    keep[i, lay$coef_idx] <- alpha
    keep[i, lay$tau_idx] <- lt
  }
  list(draws = keep, divergences = 0L, saturated = 0L)
}

#' Fit the marginal mixed-model lens
#'
#' Runs the koine-authored No-U-Turn Sampler on the random-effects *marginal*
#' model --- the random intercepts integrated out, exactly for the Gaussian LMM
#' and by adaptive Gauss-Hermite quadrature for a GLMM --- and returns an
#' [fb_fit()] of `(beta, sigma, tau)` draws. For the Poisson marginal --- whose
#' variance-component funnel a single static-metric kernel cannot traverse --- a
#' funnel-safe blocked sampler is used: `log tau` is slice-sampled, and the
#' coefficients are slice-sampled along their `tau`-resolved principal axes (a
#' Gibbs block that needs no step size and cannot stall); the Gaussian and
#' binomial marginals keep the joint NUTS kernel. The pointwise log-likelihood is
#' left `NULL` for [assemble_loglik()] to compute as the per-observation
#' integrated predictive (F3). If the pass carries no supported random intercept,
#' returns a lens-scope [abstention()] rather than a wrong fit (I1).
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
#' @param max_treedepth Maximum NUTS tree depth (the joint-NUTS Gaussian / binomial
#'   path only), capping the leapfrog budget per iteration at `2^max_treedepth`;
#'   saturation of the cap is surfaced in the provenance trace, never silently
#'   absorbed. The Poisson blocked slice kernel has no tree depth.
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
  lpg_val <- function(theta) {
    .glmm_marginal_logpost_value(theta, design, redesign, pass, lay, gh, prior,
                                 rep_pos)
  }
  # The Poisson marginal's variance-component funnel needs the blocked kernel
  # (slice on log tau, slice the coefficients along their tau-resolved principal
  # axes at fixed tau); the Gaussian and binomial marginals are well-served by the
  # joint NUTS kernel (validated, near-zero divergences).
  use_blocked <- rep_pos > 0L
  chains <- vector("list", n_chains)
  divergences <- 0L
  saturated <- 0L
  for (c in seq_len(n_chains)) {
    init <- numeric(lay$dim)
    init[lay$coef_idx] <- start + stats::rnorm(design$K, 0, 0.5)
    if (lay$has_sigma) init[lay$sigma_idx] <- stats::rnorm(1, 0, 0.5)
    init[lay$tau_idx] <- stats::rnorm(1, 0, 0.5)
    # target_accept 0.9 (above the NUTS default 0.8) and max_treedepth cap the
    # joint-NUTS Gaussian / binomial path; the Poisson blocked path needs neither.
    ch <- if (use_blocked) {
      .glmm_marginal_blocked_chain(init, lpg, lpg_val, lay, n_iter = n_iter,
                                   n_warmup = n_burn)
    } else {
      .nuts_one_chain(init, lpg, n_iter = n_iter, n_warmup = n_burn,
                      target_accept = 0.9, max_treedepth = max_treedepth)
    }
    chains[[c]] <- ch$draws
    divergences <- divergences + ch$divergences
    saturated <- saturated + ch$saturated
  }

  draws <- .glmm_marginal_to_draws(chains, design, lay, rep_pos)
  # Post-fit convergence gate (blocked Poisson path). The marginalised variance-
  # component funnel at a large random-intercept SD (tau beyond ~2, the ~8% upper
  # tail of the prior) cannot be mixed cheaply by any marginal sampler --- with the
  # random effects integrated out, the non-centred remedy is unavailable. Rather
  # than contribute a badly-mixed opinion, the lens ABSTAINS when the fit has not
  # converged: too few effective draws of any parameter (the funnel collapses
  # ess_tau to ~3 against ~100+ in the identifiable regime) or chains that disagree
  # (rhat > 1.1). This is a first-class, honest outcome (SS3) --- koine contributes
  # a calibrated Poisson opinion where the variance component is identifiable, and
  # abstains where it is not. The Gaussian/binomial joint-NUTS path mixes cleanly
  # and is not gated here (its divergences are surfaced in the trace).
  if (use_blocked) {
    conv <- posterior::summarise_draws(
      draws, ess_bulk = posterior::ess_bulk, rhat = posterior::rhat)
    ess_floor <- suppressWarnings(min(conv$ess_bulk, na.rm = TRUE))
    rhat_max <- if (n_chains >= 2L) {
      suppressWarnings(max(conv$rhat, na.rm = TRUE))
    } else {
      1
    }
    if (!is.finite(ess_floor) || ess_floor < 10 ||
          (is.finite(rhat_max) && rhat_max > 1.1)) {
      return(abstention(
        "convergence_failed", scope = "lens",
        detail = sprintf(paste0(
          "Poisson marginal under-mixed (min ess_bulk %.1f, max rhat %.2f): the ",
          "large random-intercept-SD funnel is too weakly identified to mix, so ",
          "koine abstains rather than contribute a badly-mixed opinion."),
          ess_floor, rhat_max)))
    }
  }
  # The Gaussian LMM marginal is exact; the GLMM marginal is a quadrature
  # approximation --- reported honestly per family (I1).
  exact <- if (lay$has_sigma) "exact" else "approximate"
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.1.0", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = exact,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(koine = "0.2.0"),
    trace = list(fit = if (use_blocked) {
      sprintf("glmm_marginal/blocked-slice completed (%d GH nodes)", n_nodes)
    } else {
      sprintf(paste0("glmm_marginal/NUTS completed (%d GH nodes, ",
                     "%d divergences, %d treedepth-saturated)"),
              n_nodes, divergences, saturated)
    }))
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
