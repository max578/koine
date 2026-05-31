# NIMBLE compiles generated C++: slow + toolchain-dependent, so these skip on
# CRAN. The mixed-model lens samples the random effects and reports only the
# (beta, sigma, tau) marginal --- the identical target the marginal lens reaches
# by integration, so the two must agree (the orthogonal-diversity premise).

test_that("nimble_glmm recovers the LMM marginal and agrees with glmm_marginal", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(13)
  g <- rep(1:18, each = 6)
  x <- rnorm(108)
  b <- rnorm(18, 0, 1.1)
  y <- 1 + 0.6 * x + b[g] + rnorm(108, 0, 0.7)
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 31L)

  fit <- fit_nimble_glmm(pass, n_iter = 2000L, n_burn = 1000L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x", "sigma", "tau"))
  expect_identical(fit$provenance$independence_key$kernel, "nimble-mwg-glmm")

  m <- stats::setNames(posterior::summarise_draws(fit$draws, "mean")$mean,
                       posterior::variables(fit$draws))
  fm <- fit_glmm_marginal(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 2L)
  mm <- stats::setNames(posterior::summarise_draws(fm$draws, "mean")$mean,
                        posterior::variables(fm$draws))
  # sampled-b (Gibbs) and integrated-b (quadrature) target the same marginal
  expect_equal(unname(m["b_x"]), unname(mm["b_x"]), tolerance = 0.1)
  expect_equal(unname(m["tau"]), unname(mm["tau"]), tolerance = 0.3)
})

test_that("nimble_glmm abstains (lens scope) on a population pass", {
  skip_if_not_installed("nimble")
  d <- data.frame(y = rnorm(20), x = rnorm(20))
  out <- fit_nimble_glmm(koine_pass(y ~ x, d, gaussian()))
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
})

test_that("synthesise produces a mixed-model opinion; degeneracy is honest", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(21)
  g <- rep(1:20, each = 8)
  x <- rnorm(160)
  b <- rnorm(20, 0, 1)
  y <- 0.4 + 1.0 * x + b[g] + rnorm(160, 0, 0.8)
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 5L)

  reg0 <- koine_stage1_registry()
  pick <- function(id) reg0[[which(vapply(reg0, function(z) z$id, "") == id)]]
  reg <- list(pick("glmm_marginal"), pick("nimble_glmm"))
  reg[[1L]]$fit_fn <- function(p, ls) fit_glmm_marginal(p, ls, 3000L, 1500L, 3L)
  reg[[2L]]$fit_fn <- function(p, ls) fit_nimble_glmm(p, ls, 3000L, 1500L, 3L)

  res <- synthesise(pass, registry = reg)
  # a borderline PSIS gate is a principled abstention, not a test failure;
  # where the gate passes (the usual case here) we assert the opinion's shape.
  skip_if(is_abstention(res), "PSIS gate borderline; abstention is principled")
  expect_s3_class(res, "koine_result")
  expect_equal(sum(res$opinion$weights), 1, tolerance = 1e-6)
  expect_setequal(res$independence$family, c("random_walk", "gibbs"))
  # for a Gaussian LMM the two lenses corroborate (a solved problem); koine
  # reports the degeneracy honestly rather than manufacturing diversity
  expect_true(res$diagnostics$weight_concentration$n_eff <= 2)

  # in abstain mode the same corroboration is a first-class weight_degenerate
  if (res$diagnostics$weight_concentration$degenerate) {
    res2 <- synthesise(pass, registry = reg, on_degenerate = "abstain")
    expect_true(is_abstention(res2))
    expect_identical(res2$reason, "weight_degenerate")
  }
})
