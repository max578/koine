## ----------------------------------------------------------------------------
## NUTS kernel --- No-U-Turn Sampler with dual averaging (koine-authored)
##
## A self-contained, base-R Hamiltonian Monte Carlo sampler in the
## `hamiltonian / nuts` independence family. It is the gradient-informed kernel
## the marginal mixed-model lens (`fit_glmm_marginal`) rides: the random
## intercepts are integrated out (adaptive Gauss-Hermite, `re_marginal.R`) and
## NUTS samples the smooth, low-dimensional marginal target `(beta, log tau)`.
##
## Why HMC, not the random walk: the confounded count-intercept ridge of the
## Poisson marginal mixes far too slowly under a random-walk Metropolis kernel
## (SBC-failing intercept ESS even after the level reparameterisation). A
## gradient-informed trajectory traverses the ridge in one sweep, so the marginal
## posterior is sampled to calibration where the random walk could not.
##
## The implementation follows Hoffman & Gelman (2014, JMLR; "The No-U-Turn
## Sampler"): the recursive slice-sampler tree builder (Algorithm 3) with dual-
## averaging step-size adaptation (Algorithm 6, Nesterov 2009). It adds a dense
## mass matrix estimated from a warm-up window (standard Stan practice), so any
## residual posterior correlation the reparameterisation leaves is absorbed into
## the metric rather than left for the sampler to crawl through. The kernel is
## generic over a `theta -> list(value, grad)` target, so any koine lens with an
## evaluable marginal gradient reuses it (mirroring `.dram_one_chain`).
##
## It shares no code with Stan or any HMC engine --- a genuinely independent
## Hamiltonian witness, distinct from the brms/Stan engine koine triangulates
## against and from the random-walk DRAM/DREAM lenses in the same battery.
## ----------------------------------------------------------------------------

# Maximum tolerated drop in the joint log-density along a trajectory before a
# leapfrog step is declared divergent (Hoffman & Gelman's Delta_max = 1000).
.nuts_delta_max <- 1000

# One leapfrog step of size `eps` under the dense inverse-mass `m_inv` (the
# position metric). Carries the gradient forward so the next step reuses it.
# `lpg` is the target: theta -> list(value, grad).
.nuts_leapfrog <- function(state, eps, m_inv, lpg) {
  p_half <- state$p + (eps / 2) * state$grad
  theta <- state$theta + eps * as.numeric(m_inv %*% p_half)
  vg <- lpg(theta)
  p <- p_half + (eps / 2) * vg$grad
  list(theta = theta, p = p, grad = vg$grad, logp = vg$value)
}

# Joint (negative Hamiltonian): log target + log kinetic = logp - 0.5 p' m_inv p.
# Sampling stays in this joint space so the slice variable is exact.
.nuts_joint <- function(state, m_inv) {
  state$logp - 0.5 * sum(state$p * as.numeric(m_inv %*% state$p))
}

# Generalised no-U-turn criterion under the metric: continue (TRUE) while the
# trajectory still makes forward progress at both ends, measured with the
# velocities v = m_inv p (the mass-matrix-correct inner product).
.nuts_no_uturn <- function(theta_minus, theta_plus, p_minus, p_plus, m_inv) {
  dtheta <- theta_plus - theta_minus
  v_minus <- as.numeric(m_inv %*% p_minus)
  v_plus <- as.numeric(m_inv %*% p_plus)
  sum(dtheta * v_minus) >= 0 && sum(dtheta * v_plus) >= 0
}

# Heuristic for a reasonable starting step size (Hoffman & Gelman Algorithm 4):
# from a random momentum, double or halve `eps` until a single leapfrog crosses
# acceptance probability 0.5, so dual averaging starts in the right regime.
.nuts_find_epsilon <- function(theta0, grad0, logp0, m_inv, m_chol, lpg) {
  eps <- 1
  p0 <- as.numeric(crossprod(m_chol, stats::rnorm(length(theta0))))
  state <- list(theta = theta0, p = p0, grad = grad0, logp = logp0)
  joint0 <- .nuts_joint(state, m_inv)
  prop <- .nuts_leapfrog(state, eps, m_inv, lpg)
  log_ratio <- .nuts_joint(prop, m_inv) - joint0
  a <- if (is.finite(log_ratio) && log_ratio > log(0.5)) 1 else -1
  it <- 0L
  while (is.finite(log_ratio) && a * log_ratio > -a * log(2) && it < 100L) {
    eps <- eps * 2^a
    prop <- .nuts_leapfrog(state, eps, m_inv, lpg)
    log_ratio <- .nuts_joint(prop, m_inv) - joint0
    it <- it + 1L
  }
  eps
}

# Recursive tree builder (Hoffman & Gelman Algorithm 3, slice variant). Doubles
# the trajectory in direction `dir` (+/-1) to depth `depth`, returning the two
# leaf edges, a uniformly-sampled candidate from the valid (in-slice, non-
# divergent) subset, the in-slice count `n`, the continue flag `s`, and the
# running acceptance statistic (`alpha`, `n_alpha`) dual averaging consumes.
.nuts_build_tree <- function(state, log_u, dir, depth, eps, joint0, m_inv,
                             lpg) {
  if (depth == 0L) {
    leaf <- .nuts_leapfrog(state, dir * eps, m_inv, lpg)
    joint <- .nuts_joint(leaf, m_inv)
    # a non-finite joint (extreme theta overflowing the marginal) is a hard
    # divergence: out of slice, stop, never selected. Keeps n / s in {0, 1}
    # rather than NA, so the tree-combine tests stay well-defined.
    if (!is.finite(joint)) {
      return(list(minus = leaf, plus = leaf, prop = leaf, n = 0L, s = 0L,
                  alpha = 0, n_alpha = 1L, diverged = 1L))
    }
    n <- as.integer(log_u <= joint)
    s <- as.integer(joint - log_u > -.nuts_delta_max)        # not divergent
    alpha <- min(1, exp(joint - joint0))
    if (!is.finite(alpha)) alpha <- 0
    return(list(minus = leaf, plus = leaf, prop = leaf, n = n, s = s,
                alpha = alpha, n_alpha = 1L, diverged = 1L - s))
  }
  # build the first (inner) subtree, then expand outward in the same direction
  sub <- .nuts_build_tree(state, log_u, dir, depth - 1L, eps, joint0, m_inv,
                          lpg)
  if (sub$s == 0L) {
    return(sub)                                              # already stopped
  }
  edge <- if (dir == -1L) sub$minus else sub$plus
  sub2 <- .nuts_build_tree(edge, log_u, dir, depth - 1L, eps, joint0, m_inv,
                           lpg)
  if (dir == -1L) {
    minus <- sub2$minus
    plus <- sub$plus
  } else {
    minus <- sub$minus
    plus <- sub2$plus
  }
  # progressive (memory-light) uniform sampling over the in-slice points
  prop <- sub$prop
  n_total <- sub$n + sub2$n
  if (sub2$n > 0L && stats::runif(1) < sub2$n / max(n_total, 1L)) {
    prop <- sub2$prop
  }
  s <- sub2$s *
    as.integer(.nuts_no_uturn(minus$theta, plus$theta, minus$p, plus$p, m_inv))
  list(minus = minus, plus = plus, prop = prop, n = n_total, s = s,
       alpha = sub$alpha + sub2$alpha, n_alpha = sub$n_alpha + sub2$n_alpha,
       diverged = sub$diverged + sub2$diverged)
}

# One NUTS transition from `state` at step size `eps` under metric `m_inv`
# (Cholesky `m_chol` of the mass matrix for momentum resampling). Returns the
# next state plus the trajectory's mean acceptance (for dual averaging) and a
# divergence flag.
.nuts_transition <- function(state, eps, m_inv, m_chol, max_treedepth, lpg) {
  state$p <- as.numeric(crossprod(m_chol, stats::rnorm(length(state$theta))))
  joint0 <- .nuts_joint(state, m_inv)
  if (!is.finite(joint0)) {
    # degenerate start (should not occur from a kept state); stay put.
    return(list(theta = state$theta, grad = state$grad, logp = state$logp,
                mean_accept = 0, diverged = 1L))
  }
  log_u <- joint0 - stats::rexp(1)
  minus <- state
  plus <- state
  cur <- state
  n <- 1L
  s <- 1L
  alpha <- 0
  n_alpha <- 1L
  diverged <- 0L
  for (depth in seq_len(max_treedepth)) {
    dir <- if (stats::runif(1) < 0.5) -1L else 1L
    if (dir == -1L) {
      tree <- .nuts_build_tree(minus, log_u, dir, depth - 1L, eps, joint0,
                               m_inv, lpg)
      minus <- tree$minus
    } else {
      tree <- .nuts_build_tree(plus, log_u, dir, depth - 1L, eps, joint0,
                               m_inv, lpg)
      plus <- tree$plus
    }
    alpha <- tree$alpha
    n_alpha <- tree$n_alpha
    diverged <- diverged + tree$diverged
    if (tree$s == 1L && tree$n > 0L &&
          stats::runif(1) < min(1, tree$n / n)) {
      cur <- tree$prop
    }
    n <- n + tree$n
    s <- tree$s *
      as.integer(.nuts_no_uturn(minus$theta, plus$theta, minus$p, plus$p,
                                m_inv))
    if (s == 0L) break
  }
  list(theta = cur$theta, grad = cur$grad, logp = cur$logp,
       mean_accept = alpha / n_alpha, diverged = as.integer(diverged > 0L),
       depth = depth)
}

# Dual-averaging step-size update (Nesterov 2009; Hoffman & Gelman Algorithm 6).
# Tracks the running log-step `x_bar` toward the target acceptance `delta`.
.nuts_dual_average <- function(da, accept, m, delta = 0.8, gamma = 0.05,
                               t0 = 10, kappa = 0.75) {
  da$h_bar <- (1 - 1 / (m + t0)) * da$h_bar + (1 / (m + t0)) * (delta - accept)
  log_eps <- da$mu - sqrt(m) / gamma * da$h_bar
  eta <- m^(-kappa)
  da$x_bar <- eta * log_eps + (1 - eta) * da$x_bar
  da$eps <- exp(log_eps)
  da$eps_bar <- exp(da$x_bar)
  da
}

# Run one NUTS chain over a differentiable log target `lpg` (theta ->
# list(value, grad)). Warm-up runs in two windows: the first dual-averages the
# step size under an identity metric and collects positions; a dense mass matrix
# is then estimated (regularised toward a scaled identity, Stan-style) and the
# second window re-adapts the step size under it. After warm-up both step size
# and metric are frozen and `n_iter - n_warmup` draws are kept.
#
# Returns a list with `draws` (kept x dim matrix), the frozen `eps` and inverse
# mass `m_inv`, and the post-warm-up `divergences` count (surfaced, never
# silently dropped --- a high rate flags an untrustworthy fit).
.nuts_one_chain <- function(init, lpg, n_iter, n_warmup, max_treedepth = 10L,
                            target_accept = 0.8) {
  d <- length(init)
  vg <- lpg(init)
  if (!is.finite(vg$value)) {
    # escape an impossible warm start: jitter until the target is finite
    for (try in seq_len(50L)) {
      init <- init + stats::rnorm(d, 0, 0.25)
      vg <- lpg(init)
      if (is.finite(vg$value)) break
    }
  }
  state <- list(theta = init, p = numeric(d), grad = vg$grad, logp = vg$value)

  warm1 <- n_warmup %/% 2L
  warm2 <- n_warmup - warm1
  n_keep <- n_iter - n_warmup
  keep <- matrix(NA_real_, nrow = n_keep, ncol = d)

  run_window <- function(state, m_inv, m_chol, n_win, collect) {
    eps0 <- .nuts_find_epsilon(state$theta, state$grad, state$logp, m_inv,
                               m_chol, lpg)
    da <- list(mu = log(10 * eps0), eps = eps0, eps_bar = 1, h_bar = 0,
               x_bar = 0)
    seen <- if (collect) matrix(NA_real_, nrow = n_win, ncol = d) else NULL
    for (m in seq_len(n_win)) {
      tr <- .nuts_transition(state, da$eps, m_inv, m_chol, max_treedepth, lpg)
      state <- list(theta = tr$theta, p = numeric(d), grad = tr$grad,
                    logp = tr$logp)
      da <- .nuts_dual_average(da, tr$mean_accept, m, delta = target_accept)
      if (collect) seen[m, ] <- state$theta
    }
    list(state = state, eps = da$eps_bar, seen = seen)
  }

  # window 1: identity metric, adapt step size, collect positions
  m_inv <- diag(d)
  m_chol <- diag(d)
  w1 <- run_window(state, m_inv, m_chol, warm1, collect = TRUE)
  state <- w1$state

  # estimate a dense mass matrix from window-1 positions (m_inv = posterior cov,
  # regularised toward a scaled identity so a short / collinear window stays
  # invertible), then take its Cholesky for momentum resampling.
  seen <- w1$seen[apply(w1$seen, 1L, function(r) all(is.finite(r))), ,
                  drop = FALSE]
  if (nrow(seen) > d + 1L) {
    s_cov <- stats::cov(seen)
    n_s <- nrow(seen)
    s_cov <- (n_s / (n_s + 5)) * s_cov +
      (5 / (n_s + 5)) * 1e-3 * diag(d)
    m_inv <- s_cov
  }
  m_inv <- .nuts_regularise(m_inv, d)
  m_mat <- solve(m_inv)
  m_chol <- chol(m_mat)

  # window 2: re-adapt step size under the new metric, discard draws
  w2 <- run_window(state, m_inv, m_chol, warm2, collect = FALSE)
  state <- w2$state
  eps <- w2$eps

  # sampling: frozen step size + metric
  divergences <- 0L
  saturated <- 0L                       # iterations that hit max tree depth
  for (i in seq_len(n_keep)) {
    tr <- .nuts_transition(state, eps, m_inv, m_chol, max_treedepth, lpg)
    state <- list(theta = tr$theta, p = numeric(d), grad = tr$grad,
                  logp = tr$logp)
    keep[i, ] <- state$theta
    divergences <- divergences + tr$diverged
    saturated <- saturated + as.integer(tr$depth >= max_treedepth)
  }
  list(draws = keep, eps = eps, m_inv = m_inv, divergences = divergences,
       saturated = saturated)
}

# Symmetrise and floor the eigenvalues of an estimated metric so it is a valid
# (positive-definite) inverse-mass matrix even when the warm-up covariance is
# near-singular (a tightly-identified or near-collinear marginal).
.nuts_regularise <- function(m_inv, d) {
  m_inv <- (m_inv + t(m_inv)) / 2
  e <- eigen(m_inv, symmetric = TRUE)
  vals <- pmax(e$values, 1e-8)
  e$vectors %*% (vals * t(e$vectors))
}
