test_that("gaussian family_spec matches dnorm", {
  fs <- as_family_spec(gaussian())
  y <- c(-1, 0, 2.5)
  eta <- c(0, 1, 2)
  expect_equal(
    fs$loglik_pointwise(y, eta, list(sigma = 1.3)),
    stats::dnorm(y, mean = eta, sd = 1.3, log = TRUE)
  )
  expect_identical(fs$name, "gaussian")
  expect_identical(fs$aux_params, "sigma")
})

test_that("binomial family_spec matches dbinom (Bernoulli + trials)", {
  fs <- as_family_spec(binomial())
  y <- c(0, 1, 1)
  eta <- c(-1, 0, 2)
  expect_equal(
    fs$loglik_pointwise(y, eta, list()),
    stats::dbinom(y, size = 1, prob = stats::plogis(eta), log = TRUE)
  )
  fs2 <- as_family_spec(binomial(), trials = c(5L, 5L, 5L))
  yk <- c(2, 3, 5)
  expect_equal(
    fs2$loglik_pointwise(yk, eta, list()),
    stats::dbinom(yk, size = 5, prob = stats::plogis(eta), log = TRUE)
  )
})

test_that("poisson family_spec matches dpois", {
  fs <- as_family_spec(poisson())
  y <- c(0, 2, 5)
  eta <- c(0, 0.5, 1)
  expect_equal(
    fs$loglik_pointwise(y, eta, list()),
    stats::dpois(y, lambda = exp(eta), log = TRUE)
  )
})

test_that("unsupported family aborts with a typed error", {
  expect_error(as_family_spec("weibull"), class = "koine_error_family_not_evaluable")
  expect_error(as_family_spec(Gamma()), class = "koine_error_family_not_evaluable")
})

test_that("support_check rejects out-of-support data", {
  fs <- as_family_spec(poisson())
  expect_false(fs$support_check(c(1, -2)))
  expect_false(fs$support_check(c(1, 2.5)))
  expect_true(fs$support_check(c(0, 1, 7)))
})
