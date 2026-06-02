test_that("koine_can_contribute gates the supported model corner", {
  d <- data.frame(y = rnorm(50), x = rnorm(50))
  expect_true(koine_can_contribute(y ~ x, d, gaussian()))

  # an unsupported random-effects structure -> abstention reason, not a crash
  dre <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:10, 6)))
  expect_identical(koine_can_contribute(y ~ x + (x | g), dre, gaussian()),
                   "feature_unsupported")

  # a one-lens registry cannot field a diverse opinion
  one <- list(koine_stage1_registry()[[1L]])
  expect_identical(koine_can_contribute(y ~ x, d, gaussian(), registry = one),
                   "no_diverse_battery")

  # a malformed model (column not in data) is reported, never silently fitted
  expect_identical(koine_can_contribute(y ~ nope, d, gaussian()), "pass_invalid")
})

test_that("synthesise keeps fits only on request", {
  set.seed(1)
  d <- data.frame(x = rnorm(80))
  d$y <- 1 + 0.6 * d$x + rnorm(80)
  pass <- koine_pass(y ~ x, d, gaussian(), seed = 4L)
  res_default <- synthesise(pass)
  res_keep <- synthesise(pass, keep_fits = TRUE)
  skip_if(is_abstention(res_default))
  expect_null(res_default$fits)
  expect_true(length(res_keep$fits) >= 2L)
})

test_that("koine_opinion_draws materialises the stacking mixture over shared params", {
  set.seed(1)
  d <- data.frame(x = rnorm(80))
  d$y <- 1 + 0.6 * d$x + rnorm(80)
  pass <- koine_pass(y ~ x, d, gaussian(), seed = 4L)
  res <- synthesise(pass, keep_fits = TRUE)
  skip_if(is_abstention(res))

  dr <- koine_opinion_draws(res, n_draws = 2500L)
  expect_s3_class(dr, "draws_array")
  expect_identical(posterior::ndraws(dr), 2500L)        # allocation sums exactly
  # shared params include the fixed effects every lens estimates
  expect_true(all(c("b_(Intercept)", "b_x") %in% posterior::variables(dr)))
  # the opinion's predictive-space signal rides along, never lost
  expect_true(!is.null(attr(dr, "koine_weights")))
  expect_true(is.numeric(attr(dr, "koine_elpd_stacked")))
  # the mixture posterior tracks the data's likelihood-implied estimate (the OLS
  # slope for this sample), not the generative truth --- it is a posterior, and
  # every lens agrees on it.
  b_ols <- unname(stats::coef(stats::lm(y ~ x, d))["x"])
  m <- posterior::summarise_draws(dr, "mean")
  mm <- stats::setNames(m$mean, m$variable)
  expect_equal(unname(mm["b_x"]), b_ols, tolerance = 0.1)
})

test_that("koine_opinion_draws is reproducible and passes abstentions through", {
  set.seed(1)
  d <- data.frame(x = rnorm(80))
  d$y <- 1 + 0.6 * d$x + rnorm(80)
  pass <- koine_pass(y ~ x, d, gaussian(), seed = 4L)
  res <- synthesise(pass, keep_fits = TRUE)
  skip_if(is_abstention(res))
  a <- koine_opinion_draws(res, n_draws = 1500L)
  b <- koine_opinion_draws(res, n_draws = 1500L)
  expect_identical(as.vector(a), as.vector(b))           # seeded from ir_hash

  # an abstention is returned unchanged (no opinion to represent)
  abst <- abstention("no_diverse_battery", detail = "test")
  expect_identical(koine_opinion_draws(abst), abst)
})

test_that("koine_opinion_draws blends lenses in proportion to their weights", {
  # deterministic test of the mixing mechanism: given two lenses with distinct
  # posteriors and explicit weights, the mixture mean is the weight-blended mean
  # and the draw allocation is exactly proportional. (Natural data usually
  # corroborates onto one lens or abstains on heavy tails; here we drive the
  # weights directly to exercise the blend.)
  set.seed(1)
  d <- data.frame(x = rnorm(80))
  d$y <- 1 + 0.6 * d$x + rnorm(80)
  pass <- koine_pass(y ~ x, d, gaussian(), seed = 4L)
  res <- synthesise(pass, keep_fits = TRUE)
  skip_if(is_abstention(res))
  by_id <- stats::setNames(
    lapply(res$fits, function(f) f$draws),
    vapply(res$fits, function(f) f$provenance$method, ""))
  skip_if(!all(c("dram", "lasso") %in% names(by_id)))
  mean_bx <- function(id) {
    mean(posterior::as_draws_matrix(by_id[[id]])[, "b_x"])
  }
  res$opinion$weights <- c(dram = 0.7, lasso = 0.3)       # drive the blend
  dr <- koine_opinion_draws(res, n_draws = 4000L)
  expect_identical(posterior::ndraws(dr), 4000L)
  blended <- 0.7 * mean_bx("dram") + 0.3 * mean_bx("lasso")
  got <- mean(posterior::as_draws_matrix(dr)[, "b_x"])
  expect_equal(got, blended, tolerance = 0.03)            # weight-blended mean
})

test_that("koine_opinion_draws errors clearly without retained fits", {
  set.seed(1)
  d <- data.frame(x = rnorm(60))
  d$y <- d$x + rnorm(60)
  res <- synthesise(koine_pass(y ~ x, d, gaussian(), seed = 4L))
  skip_if(is_abstention(res))
  expect_error(koine_opinion_draws(res), "keep_fits = TRUE")
})
