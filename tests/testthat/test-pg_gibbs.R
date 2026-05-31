make_binom <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- 0.5 + 1.2 * x
  y <- stats::rbinom(n, 1, stats::plogis(eta))
  data.frame(y = y, x = x)
}

test_that("pg_gibbs recovers logistic-regression coefficients", {
  skip_if_not_installed("BayesLogit")
  d <- make_binom(seed = 2L)
  pass <- koine_pass(y ~ x, d, binomial(), seed = 7L)
  fit <- fit_pg_gibbs(pass, n_iter = 2000L, n_burn = 1000L, n_chains = 2L)

  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws), c("b_(Intercept)", "b_x"))
  expect_identical(fit$provenance$independence_key$kernel, "pg-exact")
  expect_identical(fit$provenance$maturity_tier, "koine-implemented-validated")

  s <- posterior::summarise_draws(fit$draws, "mean", rhat = posterior::rhat)
  m <- stats::setNames(s$mean, s$variable)
  glm_co <- stats::coef(stats::glm(y ~ x, binomial(), data = d))
  expect_equal(unname(m["b_(Intercept)"]), unname(glm_co[1L]), tolerance = 0.15)
  expect_equal(unname(m["b_x"]), unname(glm_co[2L]), tolerance = 0.15)
  expect_true(all(s$rhat < 1.1))
})

test_that("pg_gibbs is deterministic given the pass seed", {
  skip_if_not_installed("BayesLogit")
  pass <- koine_pass(y ~ x, make_binom(seed = 3L), binomial(), seed = 9L)
  f1 <- fit_pg_gibbs(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
  f2 <- fit_pg_gibbs(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
  expect_equal(posterior::as_draws_matrix(f1$draws),
               posterior::as_draws_matrix(f2$draws))
})

test_that("pg_gibbs abstains (lens scope) on a non-binomial family", {
  skip_if_not_installed("BayesLogit")
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(20), x = rnorm(20)),
                     gaussian())
  out <- fit_pg_gibbs(pass)
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
})

test_that("pg_gibbs does not resolve for a Gaussian pass (registry gate)", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(20), x = rnorm(20)),
                     gaussian())
  res <- resolve_battery(pass)
  ids <- vapply(res$enrolled, function(L) L$id, character(1))
  expect_false("pg_gibbs" %in% ids)
})
