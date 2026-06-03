test_that("slice sampler recovers a Gaussian target", {
  set.seed(1)
  lf <- function(x) stats::dnorm(x, 2, 0.7, log = TRUE)
  x <- 0
  draws <- numeric(4000L)
  for (i in seq_along(draws)) {
    x <- koine:::.slice_1d(lf, x, w = 1)
    draws[i] <- x
  }
  expect_equal(mean(draws), 2, tolerance = 0.05)
  expect_equal(stats::sd(draws), 0.7, tolerance = 0.05)
})

test_that("slice sampler recovers a skewed (log-scale) target", {
  # x = log(g), g ~ Gamma(3, 2): the density of x is dgamma(exp(x)) * exp(x).
  set.seed(2)
  lf <- function(x) stats::dgamma(exp(x), 3, 2, log = TRUE) + x
  x <- 0
  draws <- numeric(5000L)
  for (i in seq_along(draws)) {
    x <- koine:::.slice_1d(lf, x, w = 1)
    draws[i] <- x
  }
  expect_equal(mean(exp(draws)), 1.5, tolerance = 0.05)        # Gamma(3,2) mean
})

test_that("slice sampler respects hard support bounds", {
  set.seed(3)
  lf <- function(x) stats::dnorm(x, 0, 1, log = TRUE)          # mass below 0 too
  x <- 0.5
  draws <- numeric(2000L)
  for (i in seq_along(draws)) {
    x <- koine:::.slice_1d(lf, x, w = 1, lower = 0)            # truncate at 0
    draws[i] <- x
  }
  expect_true(all(draws >= 0))
  expect_gt(mean(draws), 0.6)                                  # half-normal mean
})

test_that("slice sampler stays put on a non-finite start (I2-safe)", {
  lf <- function(x) if (x < 0) -Inf else stats::dnorm(x, 0, 1, log = TRUE)
  expect_identical(koine:::.slice_1d(lf, -1, w = 1), -1)       # degenerate, held
})

test_that("slice sampler leaves a bounded mixture invariant", {
  # a two-component target; the sampler must spend the right mass in each well.
  set.seed(4)
  lf <- function(x) {
    log(0.3 * stats::dnorm(x, -2, 0.5) + 0.7 * stats::dnorm(x, 2, 0.5))
  }
  x <- 0
  draws <- numeric(6000L)
  for (i in seq_along(draws)) {
    x <- koine:::.slice_1d(lf, x, w = 2)
    draws[i] <- x
  }
  expect_equal(mean(draws > 0), 0.7, tolerance = 0.05)         # right-well mass
})
