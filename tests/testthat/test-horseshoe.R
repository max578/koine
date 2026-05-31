sparse_pass <- function(n = 80L, p = 4L, sig = 1.5, seed = 3L, ps = 4L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  bt <- c(sig, rep(0, p - 1L))                       # one signal, rest noise
  y <- 1 + as.numeric(X %*% bt) + rnorm(n, sd = 1)
  d <- data.frame(y = y, X)
  form <- stats::as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
  koine_pass(form, d, gaussian(), seed = ps)
}

test_that("horseshoe recovers the signal and shrinks the noise coefficients", {
  skip_if_not_installed("bayesreg")
  pass <- sparse_pass(n = 100L, p = 4L, sig = 1.6)
  fit <- fit_horseshoe(pass, n_iter = 1500L, n_burn = 1000L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x1", "b_x2", "b_x3", "b_x4", "sigma"))
  expect_silent(validate_provenance(fit$provenance))

  m <- stats::setNames(posterior::summarise_draws(fit$draws, "mean")$mean,
                       posterior::summarise_draws(fit$draws, "mean")$variable)
  expect_equal(unname(m["b_x1"]), 1.6, tolerance = 0.25)   # signal recovered
  # the three noise coefficients are shrunk toward zero
  expect_true(all(abs(m[c("b_x2", "b_x3", "b_x4")]) < 0.3))
})

test_that("horseshoe abstains (lens scope) on non-Gaussian / too-few predictors", {
  skip_if_not_installed("bayesreg")
  pb <- koine_pass(y ~ x, data.frame(y = rbinom(20, 1, 0.5), x = rnorm(20)),
                   binomial())
  expect_identical(fit_horseshoe(pb)$reason, "feature_unsupported")

  p0 <- koine_pass(y ~ 1, data.frame(y = rnorm(20)), gaussian())
  out <- fit_horseshoe(p0)
  expect_true(is_abstention(out))                          # intercept-only
  expect_identical(out$reason, "feature_unsupported")
})

test_that("horseshoe vs flat-normal realises sparsity-structure diversity", {
  skip_on_cran()
  skip_if_not_installed("bayesreg")
  pass <- sparse_pass(n = 120L, p = 3L, sig = 0.6, seed = 3L, ps = 4L)
  reg0 <- koine_stage1_registry()
  pick <- function(id) reg0[[which(vapply(reg0, function(z) z$id, "") == id)]]
  two <- list(pick("dram"), pick("horseshoe"))
  two[[1L]]$fit_fn <- function(p, ls) fit_dram(p, ls, 6000L, 3000L, 2L)
  two[[2L]]$fit_fn <- function(p, ls) fit_horseshoe(p, ls, 1500L, 1000L, 2L)

  res <- synthesise(pass, registry = two)
  expect_s3_class(res, "koine_result")
  wc <- res$diagnostics$weight_concentration
  expect_gt(wc$n_eff, 1.3)                  # weights do not collapse
  expect_false(wc$degenerate)
  expect_true(all(res$opinion$weights > 0.1))
})
