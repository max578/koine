# NIMBLE compiles generated C++: these tests are slow and need a toolchain, so
# they skip on CRAN. One compile total (the recovery test); the abstention test
# returns before any model is built.

test_that("nimble_gibbs recovers the Gaussian posterior and agrees with DRAM", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(13)
  n <- 80L
  x <- rnorm(n)
  y <- 1.5 + 0.8 * x + rnorm(n, sd = 0.7)
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 31L)

  fit <- fit_nimble(pass, n_iter = 1500L, n_burn = 750L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x", "sigma"))
  expect_silent(validate_provenance(fit$provenance))
  expect_identical(fit$provenance$independence_key$family, "gibbs")

  s <- posterior::summarise_draws(fit$draws, "mean")
  m <- stats::setNames(s$mean, s$variable)
  ols <- stats::coef(stats::lm(y ~ x))
  expect_equal(unname(m["b_(Intercept)"]), unname(ols[1L]), tolerance = 0.1)
  expect_equal(unname(m["b_x"]), unname(ols[2L]), tolerance = 0.1)

  # cross-family agreement: NIMBLE-Gibbs and DRAM target the SAME posterior, so
  # their posterior means must agree (the orthogonal-diversity premise --- they
  # differ in inference family, not in target).
  fit_d <- fit_dram(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 2L)
  sd <- posterior::summarise_draws(fit_d$draws, "mean")
  md <- stats::setNames(sd$mean, sd$variable)
  expect_equal(unname(m["b_x"]), unname(md["b_x"]), tolerance = 0.08)
})

test_that("nimble_gibbs abstains (lens scope) on an unsupported link", {
  skip_if_not_installed("nimble")
  d <- data.frame(y = rpois(20, 3), x = rnorm(20))
  pass <- koine_pass(y ~ x, d, gaussian(link = "log"))
  out <- fit_nimble(pass)
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
  expect_identical(out$scope, "lens")
})

test_that("synthesise yields a multi-family stacked opinion (Gaussian battery)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if_not_installed("BayesianTools")
  set.seed(21)
  n <- 70L
  x <- rnorm(n)
  y <- 0.4 + 1.1 * x + rnorm(n, sd = 0.8)
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 5L)

  # an explicit four-lens registry (small-iter wrappers), decoupled from the
  # full battery so the test survives further lenses being added
  reg0 <- koine_stage1_registry()
  pick <- function(id) reg0[[which(vapply(reg0, function(z) z$id, "") == id)]]
  reg <- list(pick("dram"), pick("nimble_gibbs"), pick("dream_zs"),
              pick("robust_t"))
  reg[[1L]]$fit_fn <- function(p, ls) fit_dram(p, ls, 1500L, 750L, 2L)
  reg[[2L]]$fit_fn <- function(p, ls) fit_nimble(p, ls, 1500L, 750L, 2L)
  reg[[3L]]$fit_fn <- function(p, ls) fit_dream(p, ls, 8000L)
  reg[[4L]]$fit_fn <- function(p, ls) fit_robust_t(p, ls, nu = 4,
                                                   n_iter = 1500L,
                                                   n_burn = 750L, n_chains = 2L)

  res <- synthesise(pass, registry = reg)
  expect_s3_class(res, "koine_result")
  expect_setequal(names(res$opinion$weights),
                  c("dram", "nimble_gibbs", "dream_zs", "robust_t"))
  expect_equal(sum(res$opinion$weights), 1, tolerance = 1e-6)
  # three distinct coarse families (robust_t shares random_walk with dram)
  expect_setequal(res$independence$family,
                  c("random_walk", "gibbs", "population"))
})
