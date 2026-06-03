test_that("glmm_marginal recovers a Gaussian LMM posterior", {
  set.seed(42)
  g <- rep(1:20, each = 6)
  x <- rnorm(120)
  b <- rnorm(20, 0, 1.1)
  y <- 1 + 0.5 * x + b[g] + rnorm(120, 0, 0.7)
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 7L)

  fit <- fit_glmm_marginal(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 3L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$exactness, "exact")        # Gaussian marginal is exact
  s <- posterior::summarise_draws(fit$draws, "mean", rhat = posterior::rhat)
  m <- stats::setNames(s$mean, s$variable)
  expect_equal(unname(m["b_x"]), 0.5, tolerance = 0.15)
  expect_true(m[["tau"]] > 0.7 && m[["tau"]] < 1.8)
  expect_true(m[["sigma"]] > 0.5 && m[["sigma"]] < 0.9)
  expect_true(all(s$rhat < 1.1))
})

test_that("glmm_marginal draws carry b_<col> / sigma / tau and complete provenance", {
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:10, 6)))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 3L)
  fit <- fit_glmm_marginal(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x", "sigma", "tau"))
  expect_silent(validate_provenance(fit$provenance))
  expect_identical(fit$provenance$independence_key$kernel,
                   "glmm-marginal-quadrature")
  # the lens rides the koine-authored Hamiltonian kernel (distinct inference
  # family from the random-walk DRAM/DREAM lenses in the same battery).
  expect_identical(fit$provenance$independence_key$family, "hamiltonian")
})

test_that("glmm_marginal on a GLMM is reported as an approximate (quadrature) fit", {
  set.seed(11)
  g <- rep(1:25, each = 8)
  x <- rnorm(200)
  b <- rnorm(25, 0, 1)
  y <- rbinom(200, 1, plogis(0.3 + 0.8 * x + b[g]))
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, binomial(), seed = 7L)
  fit <- fit_glmm_marginal(pass, n_iter = 2000L, n_burn = 1000L, n_chains = 2L)
  expect_identical(fit$exactness, "approximate")
  m <- posterior::summarise_draws(fit$draws, "mean")
  bx <- m$mean[m$variable == "b_x"]
  expect_equal(bx, 0.8, tolerance = 0.3)
})

test_that("glmm_marginal fits a Poisson GLMM with a well-mixed intercept", {
  set.seed(5)
  g <- rep(1:20, each = 8)
  x <- rnorm(160)
  b <- rnorm(20, 0, 0.7)
  y <- rpois(160, exp(0.3 + 0.4 * x + b[g]))
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, poisson(), seed = 7L)
  fit <- fit_glmm_marginal(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 3L)
  expect_identical(fit$status, "contributed")
  s <- posterior::summarise_draws(fit$draws, "mean", rhat = posterior::rhat,
                                  ess_bulk = posterior::ess_bulk)
  m <- stats::setNames(s$mean, s$variable)
  expect_equal(unname(m["b_x"]), 0.4, tolerance = 0.25)
  expect_true(all(s$rhat < 1.1))
  # the Hamiltonian kernel traverses the confounded count-intercept ridge in one
  # trajectory: the intercept mixes near-freely (the random-walk kernel stalled
  # it at ~5-12% efficiency, the SBC-failing regime).
  expect_gt(s$ess_bulk[s$variable == "b_(Intercept)"], 800)
})

test_that("the Poisson marginal rides the blocked slice-within-NUTS kernel", {
  d <- data.frame(y = rpois(60, 2), x = rnorm(60), g = factor(rep(1:10, 6)))
  pass <- koine_pass(y ~ x + (1 | g), d, poisson(), seed = 1L)
  fit <- fit_glmm_marginal(pass, n_iter = 400L, n_burn = 200L, n_chains = 1L)
  # Poisson uses the funnel-safe blocked kernel; the Gaussian/binomial marginals
  # keep the joint NUTS (asserted by the trace tag).
  expect_match(fit$provenance$trace$fit, "blocked-slice")
  dpass <- koine_pass(y ~ x + (1 | g),
                      data.frame(y = rnorm(60), x = rnorm(60),
                                 g = factor(rep(1:10, 6))),
                      gaussian(), seed = 1L)
  dfit <- fit_glmm_marginal(dpass, n_iter = 400L, n_burn = 200L, n_chains = 1L)
  expect_match(dfit$provenance$trace$fit, "glmm_marginal/NUTS")
  expect_no_match(dfit$provenance$trace$fit, "blocked")
})

test_that("the value-only marginal log-posterior matches the value+gradient path", {
  # the blocked slice kernel rides a value-only marginal (skipping the score-
  # identity gradient quadrature); it must equal the value the gradient path
  # returns, or the slice would target a different density.
  for (fam in c("gaussian", "binomial", "poisson")) {
    set.seed(match(fam, c("gaussian", "binomial", "poisson")))
    g <- rep(1:12, each = 6)
    x <- rnorm(72)
    b <- rnorm(12, 0, 0.9)
    y <- switch(fam,
                gaussian = 0.4 + 0.5 * x + b[g] + rnorm(72, 0, 0.8),
                binomial = rbinom(72, 1, plogis(0.2 + 0.6 * x + b[g])),
                poisson  = rpois(72, exp(0.2 + 0.4 * x + b[g])))
    d <- data.frame(y = y, x = x, g = factor(g))
    pass <- koine_pass(y ~ x + (1 | g), d,
                       get(fam, mode = "function")(), seed = 5L)
    design <- koine:::.pass_design(pass)
    redesign <- koine:::.pass_re_design(pass)
    lay <- koine:::.glmm_marginal_layout(pass, design$K)
    gh <- koine:::.gh_rule(20L)
    prior <- koine:::.re_prior(fam)
    rep_pos <- if (fam == "poisson") match("(Intercept)", colnames(design$X)) else 0L
    theta <- numeric(lay$dim)
    theta[lay$coef_idx] <- c(0.3, 0.5)
    if (lay$has_sigma) theta[lay$sigma_idx] <- log(0.8)
    theta[lay$tau_idx] <- log(0.9)
    v_grad <- koine:::.glmm_marginal_logpost_grad(
      theta, design, redesign, pass, lay, gh, prior, rep_pos)$value
    v_only <- koine:::.glmm_marginal_logpost_value(
      theta, design, redesign, pass, lay, gh, prior, rep_pos)
    expect_equal(v_only, v_grad, tolerance = 1e-9)
  }
})

test_that("the tau-adaptive knot axes span a positive-definite covariance", {
  set.seed(8)
  g <- rep(1:15, each = 6)
  x <- rnorm(90)
  b <- rnorm(15, 0, 1.5)
  y <- rpois(90, exp(0.3 + 0.4 * x + b[g]))
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, poisson(), seed = 3L)
  design <- koine:::.pass_design(pass)
  redesign <- koine:::.pass_re_design(pass)
  lay <- koine:::.glmm_marginal_layout(pass, design$K)
  gh <- koine:::.gh_rule(20L)
  prior <- koine:::.re_prior("poisson")
  rep_pos <- match("(Intercept)", colnames(design$X))
  lpg <- function(th) koine:::.glmm_marginal_logpost_grad(
    th, design, redesign, pass, lay, gh, prior, rep_pos)
  knots <- koine:::.glmm_build_tau_knots(lpg, lay, c(log(mean(y)), 0))
  # the scaled principal axes (columns of `axes`) span a symmetric positive-
  # definite conditional covariance (axes axes' = cov), tau-resolved: a larger tau
  # gives the weakly-identified intercept a larger conditional variance.
  cov_of <- function(ax) tcrossprod(ax)
  for (ax in knots$axes) {
    cv <- cov_of(ax)
    expect_true(isSymmetric(cv, tol = 1e-6))
    expect_true(all(eigen(cv, symmetric = TRUE)$values > 0))
  }
  v_small <- cov_of(knots$axes[[1L]])[1L, 1L]              # small-tau knot
  v_large <- cov_of(knots$axes[[length(knots$axes)]])[1L, 1L]
  expect_gt(v_large, v_small)                             # intercept variance grows in tau
})

test_that("the blocked kernel keeps a funnel-prone Poisson fit calibrated", {
  # a moderately large random-intercept SD (the funnel-prone regime that broke the
  # joint static-metric kernel: ~50% divergences, ess_tau ~ 3, and an independence
  # proposal froze at zero acceptance). The double-slice blocked kernel must mix
  # the variance component (a non-trivial effective sample size) and recover the
  # coefficient. Modest iterations keep the guard affordable.
  set.seed(2026)
  g <- rep(1:15, each = 6)
  x <- rnorm(90)
  b <- rnorm(15, 0, 1.4)
  y <- rpois(90, exp(0.3 + 0.5 * x + b[g]))
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, poisson(), seed = 11L)
  fit <- fit_glmm_marginal(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
  expect_identical(fit$status, "contributed")
  s <- posterior::summarise_draws(fit$draws, "mean", ess_bulk = posterior::ess_bulk)
  m <- stats::setNames(s$mean, s$variable)
  ess_tau <- s$ess_bulk[s$variable == "tau"]
  expect_gt(ess_tau, 20)                            # was ~3 for the failing kernels
  expect_equal(unname(m["b_x"]), 0.5, tolerance = 0.35)
})

test_that("glmm_marginal abstains on a population pass (lens scope)", {
  d <- data.frame(y = rnorm(30), x = rnorm(30))
  fit <- fit_glmm_marginal(koine_pass(y ~ x, d, gaussian()))
  expect_true(is_abstention(fit))
  expect_identical(fit$reason, "feature_unsupported")
})

test_that("the integrated predictive keeps the PSIS gate reliable for an LMM", {
  set.seed(1)
  g <- rep(1:15, each = 6)
  x <- rnorm(90)
  b <- rnorm(15, 0, 1)
  y <- 0.5 * x + b[g] + rnorm(90, 0, 0.7)
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 4L)
  fit <- fit_glmm_marginal(pass, n_iter = 1500L, n_burn = 750L, n_chains = 2L)
  ll <- assemble_loglik(list(fit), pass)
  expect_equal(ncol(ll[[1L]]$ll), nrow(d))         # per-observation predictive
  gate <- gate_psis(ll)
  expect_length(gate$dropped, 0L)                  # k below the threshold
})
