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
