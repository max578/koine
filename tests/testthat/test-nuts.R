test_that("leapfrog is volume-preserving and reversible", {
  # a quadratic target: grad = -A theta. One forward leapfrog then one with
  # negated momentum must return to the start (time-reversibility).
  a_mat <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  lpg <- function(th) list(value = -0.5 * sum(th * (a_mat %*% th)),
                           grad = -as.numeric(a_mat %*% th))
  m_inv <- diag(2)
  th0 <- c(0.3, -0.7)
  vg <- lpg(th0)
  s0 <- list(theta = th0, p = c(0.4, 0.2), grad = vg$grad, logp = vg$value)
  fwd <- koine:::.nuts_leapfrog(s0, 0.1, m_inv, lpg)
  back <- koine:::.nuts_leapfrog(
    list(theta = fwd$theta, p = -fwd$p, grad = fwd$grad, logp = fwd$logp),
    0.1, m_inv, lpg)
  expect_equal(back$theta, th0, tolerance = 1e-10)
  expect_equal(-back$p, s0$p, tolerance = 1e-10)
})

test_that("NUTS recovers a correlated bivariate Gaussian", {
  set.seed(1)
  sig <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
  prec <- solve(sig)
  lpg <- function(th) list(value = -0.5 * sum(th * (prec %*% th)),
                           grad = -as.numeric(prec %*% th))
  ch <- koine:::.nuts_one_chain(c(0, 0), lpg, n_iter = 4000L, n_warmup = 1500L)
  d <- ch$draws
  expect_equal(colMeans(d), c(0, 0), tolerance = 0.1)
  expect_equal(apply(d, 2L, stats::sd), c(1, 1), tolerance = 0.1)
  expect_equal(stats::cor(d[, 1], d[, 2]), 0.9, tolerance = 0.05)
  expect_lt(ch$divergences, 25L)            # a smooth Gaussian rarely diverges
})

test_that("the dense mass matrix absorbs wildly different scales", {
  set.seed(2)
  sds <- c(0.02, 1, 50)
  lpg <- function(th) list(value = -0.5 * sum((th / sds)^2),
                           grad = -th / sds^2)
  d <- koine:::.nuts_one_chain(c(0, 0, 0), lpg, n_iter = 4000L,
                               n_warmup = 1500L)$draws
  expect_equal(apply(d, 2L, stats::sd) / sds, c(1, 1, 1), tolerance = 0.15)
})

test_that("NUTS is reproducible under a fixed seed", {
  lpg <- function(th) list(value = -0.5 * sum(th^2), grad = -th)
  set.seed(99)
  a <- koine:::.nuts_one_chain(c(0.5, -0.5), lpg, n_iter = 800L,
                               n_warmup = 400L)$draws
  set.seed(99)
  b <- koine:::.nuts_one_chain(c(0.5, -0.5), lpg, n_iter = 800L,
                               n_warmup = 400L)$draws
  expect_identical(a, b)
})

test_that("a non-finite target value is handled as a divergence, not an error", {
  # half the space is forbidden (value -Inf); the sampler must stay in the
  # allowed half and never error.
  lpg <- function(th) {
    if (th[1L] < 0) {
      return(list(value = -Inf, grad = c(0, 0)))
    }
    list(value = -0.5 * sum(th^2), grad = -th)
  }
  set.seed(3)
  d <- koine:::.nuts_one_chain(c(1, 0), lpg, n_iter = 1500L,
                               n_warmup = 750L)$draws
  expect_true(all(is.finite(d)))
  expect_true(all(d[, 1L] >= 0))
})
