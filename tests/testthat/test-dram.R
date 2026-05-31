test_that("DRAM recovers a Gaussian linear-model posterior", {
  set.seed(7)
  n <- 80L
  x <- rnorm(n)
  y <- 1.5 + 0.8 * x + rnorm(n, sd = 0.7)
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 11L)

  fit <- fit_dram(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 3L)
  expect_s3_class(fit, "koine_fb_fit")

  s <- posterior::summarise_draws(fit$draws, "mean", rhat = posterior::rhat)
  m <- stats::setNames(s$mean, s$variable)
  ols <- stats::coef(stats::lm(y ~ x))

  expect_equal(unname(m["b_(Intercept)"]), unname(ols[1L]), tolerance = 0.1)
  expect_equal(unname(m["b_x"]), unname(ols[2L]), tolerance = 0.1)
  expect_true(m[["sigma"]] > 0.4 && m[["sigma"]] < 1.0)
  expect_true(all(s$rhat < 1.1))           # chains mixed
})

test_that("DRAM is deterministic given the pass seed", {
  pass <- koine_pass(y ~ x,
                     data.frame(y = rnorm(40), x = rnorm(40)),
                     gaussian(), seed = 99L)
  f1 <- fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
  f2 <- fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
  expect_equal(posterior::as_draws_matrix(f1$draws),
               posterior::as_draws_matrix(f2$draws))
})

test_that("DRAM draws carry the expected variable names + complete provenance", {
  pass <- koine_pass(y ~ x,
                     data.frame(y = rnorm(30), x = rnorm(30)),
                     gaussian(), seed = 3L)
  fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
  vars <- posterior::variables(fit$draws)
  expect_setequal(vars, c("b_(Intercept)", "b_x", "sigma"))
  expect_silent(validate_provenance(fit$provenance))
  expect_identical(fit$provenance$ir_hash, pass$ir_hash)
})
