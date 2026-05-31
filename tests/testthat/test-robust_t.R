test_that("robust_t recovers coefficients and is outlier-robust", {
  set.seed(7)
  n <- 100L
  x <- rnorm(n)
  y <- 1 + 0.8 * x + rnorm(n, sd = 0.5)
  y[c(1L, 2L, 3L)] <- y[c(1L, 2L, 3L)] + c(10, -10, 10)   # gross outliers
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 4L)

  fit <- fit_robust_t(pass, n_iter = 2000L, n_burn = 1000L, n_chains = 2L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "b_x", "sigma"))
  expect_silent(validate_provenance(fit$provenance))

  # native Student-t log-likelihood stored as array[chain, draw, obs]
  expect_false(is.null(fit$log_likelihood))
  expect_equal(dim(fit$log_likelihood)[3L], n)

  # robust to the outliers: slope still recovered (a Gaussian fit would not)
  m <- stats::setNames(
    posterior::summarise_draws(fit$draws, "mean")$mean,
    posterior::summarise_draws(fit$draws, "mean")$variable)
  expect_equal(unname(m["b_x"]), 0.8, tolerance = 0.2)
})

test_that("robust_t native log-lik is Student-t, and flows through assemble", {
  set.seed(3)
  n <- 40L
  x <- rnorm(n)
  pass <- koine_pass(y ~ x, data.frame(y = 1 + 0.8 * x + rnorm(n), x = x),
                     gaussian(), seed = 5L)
  fit <- fit_robust_t(pass, nu = 4, n_iter = 800L, n_burn = 400L,
                      n_chains = 2L)
  ll <- assemble_loglik(list(fit), pass)
  expect_named(ll, "robust_t")
  m <- ll[["robust_t"]]$ll
  expect_equal(ncol(m), n)
  expect_identical(colnames(m), pass$obs_dim_names$obs_id)

  # the assembled value must equal the Student-t density, NOT a Gaussian one
  dm <- unclass(posterior::as_draws_matrix(fit$draws))
  X <- model.matrix(y ~ x, pass$data)
  beta1 <- as.numeric(dm[1L, c("b_(Intercept)", "b_x")])
  mu1 <- as.numeric(X %*% beta1)
  s1 <- dm[1L, "sigma"]
  expect_equal(
    unname(m[1L, 1L]),
    unname(stats::dt((pass$data$y[1L] - mu1[1L]) / s1, df = 4, log = TRUE) -
             log(s1)))
})

test_that("robust_t abstains (lens scope) on a non-Gaussian response", {
  pass <- koine_pass(y ~ x,
                     data.frame(y = rbinom(20, 1, 0.5), x = rnorm(20)),
                     binomial())
  out <- fit_robust_t(pass)
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
})
