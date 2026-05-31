# DREAM (BayesianTools DREAMzs) --- the `dream_zs` lens. No C++ compile, but the
# ensemble sampler is moderately heavy, so the recovery test skips on CRAN.

test_that("dream_zs recovers the Gaussian linear-model posterior", {
  skip_on_cran()
  skip_if_not_installed("BayesianTools")
  set.seed(4)
  n <- 70L
  x <- rnorm(n)
  y <- 1.5 + 0.8 * x + rnorm(n, sd = 0.7)
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 17L)

  fit <- fit_dream(pass, n_iter = 12000L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x", "sigma"))
  expect_identical(fit$provenance$independence_key$family, "population")
  expect_identical(fit$exactness, "approximate")

  s <- posterior::summarise_draws(fit$draws, "mean")
  m <- stats::setNames(s$mean, s$variable)
  ols <- stats::coef(stats::lm(y ~ x))
  expect_equal(unname(m["b_(Intercept)"]), unname(ols[1L]), tolerance = 0.12)
  expect_equal(unname(m["b_x"]), unname(ols[2L]), tolerance = 0.12)
  expect_true(m[["sigma"]] > 0.4 && m[["sigma"]] < 1.1)
})

test_that("dream_zs abstains (lens scope) on an unsupported link", {
  skip_if_not_installed("BayesianTools")
  pass <- koine_pass(y ~ x, data.frame(y = rpois(20, 3), x = rnorm(20)),
                     poisson(link = "sqrt"))
  out <- fit_dream(pass)
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
})
