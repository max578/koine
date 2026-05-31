shrink_pass <- function(n = 100L, p = 4L, sig = 1.5, seed = 3L, ps = 4L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- 1 + sig * X[, 1L] + rnorm(n, sd = 1)          # one signal, rest noise
  d <- data.frame(y = y, X)
  form <- stats::as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
  koine_pass(form, d, gaussian(), seed = ps)
}

test_that("lasso recovers the signal and shrinks the noise", {
  skip_if_not_installed("bayesreg")
  pass <- shrink_pass(sig = 1.6)
  fit <- fit_lasso(pass, n_iter = 1500L, n_burn = 1000L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$provenance$independence_key$kernel, "bayesreg-lasso")
  m <- stats::setNames(posterior::summarise_draws(fit$draws, "mean")$mean,
                       posterior::summarise_draws(fit$draws, "mean")$variable)
  expect_equal(unname(m["b_x1"]), 1.6, tolerance = 0.25)
  expect_true(all(abs(m[c("b_x2", "b_x3", "b_x4")]) < 0.35))
})

test_that("ridge recovers the signal (mild shrinkage)", {
  skip_if_not_installed("bayesreg")
  pass <- shrink_pass(sig = 1.6)
  fit <- fit_ridge(pass, n_iter = 1500L, n_burn = 1000L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$provenance$independence_key$kernel, "bayesreg-ridge")
  m <- stats::setNames(posterior::summarise_draws(fit$draws, "mean")$mean,
                       posterior::summarise_draws(fit$draws, "mean")$variable)
  expect_equal(unname(m["b_x1"]), 1.6, tolerance = 0.25)
})

test_that("shrinkage lenses abstain (lens scope) on a non-Gaussian response", {
  skip_if_not_installed("bayesreg")
  pb <- koine_pass(y ~ x, data.frame(y = rbinom(20, 1, 0.5), x = rnorm(20)),
                   binomial())
  expect_identical(fit_lasso(pb)$reason, "feature_unsupported")
  expect_identical(fit_ridge(pb)$reason, "feature_unsupported")
})

test_that("the shrinkage spectrum (flat + lasso + horseshoe) realises diversity", {
  skip_on_cran()
  skip_if_not_installed("bayesreg")
  pass <- shrink_pass(n = 120L, p = 3L, sig = 0.6, seed = 3L, ps = 4L)
  reg0 <- koine_stage1_registry()
  pick <- function(id) reg0[[which(vapply(reg0, function(z) z$id, "") == id)]]
  spectrum <- list(pick("dram"), pick("lasso"), pick("horseshoe"))
  spectrum[[1L]]$fit_fn <- function(p, ls) fit_dram(p, ls, 6000L, 3000L, 2L)
  spectrum[[2L]]$fit_fn <- function(p, ls) fit_lasso(p, ls, 1500L, 1000L, 2L)
  spectrum[[3L]]$fit_fn <- function(p, ls) fit_horseshoe(p, ls, 1500L, 1000L, 2L)

  res <- synthesise(pass, registry = spectrum)
  expect_s3_class(res, "koine_result")
  expect_equal(res$opinion$n_components, 3L)
  expect_gt(res$diagnostics$weight_concentration$n_eff, 1.3)
})
