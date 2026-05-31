test_that("Gauss-Hermite rule integrates polynomials against exp(-x^2)", {
  gh <- koine:::.gh_rule(20L)
  expect_equal(sum(gh$weights), sqrt(pi))
  expect_equal(sum(gh$weights * gh$nodes^2), sqrt(pi) / 2)
  expect_equal(sum(gh$weights * gh$nodes^4), 3 * sqrt(pi) / 4)
})

test_that("Gaussian integrated predictive matches direct integration", {
  fs <- as_family_spec(gaussian())
  sigma <- 0.7
  tau <- 1.3
  ana <- koine:::.re_integrated_loglik(
    fs, y = 0.4, eta_fixed = 0.1, tau = tau, aux = list(sigma = sigma))
  brute <- log(stats::integrate(
    function(b) stats::dnorm(0.4, 0.1 + b, sigma) * stats::dnorm(b, 0, tau),
    -Inf, Inf)$value)
  expect_equal(ana, brute, tolerance = 1e-6)
})

test_that("Gaussian integrated predictive equals N(y, xbeta, sqrt(sigma^2+tau^2))", {
  fs <- as_family_spec(gaussian())
  y <- c(-1, 0.5, 2)
  eta <- c(0.2, -0.3, 1)
  got <- koine:::.re_integrated_loglik(fs, y, eta, tau = 0.9,
                                                 aux = list(sigma = 0.6))
  want <- stats::dnorm(y, eta, sqrt(0.6^2 + 0.9^2), log = TRUE)
  expect_equal(got, want)
})

test_that("binomial / poisson integrated predictive match quadrature integration", {
  gh <- koine:::.gh_rule(24L)
  fb <- as_family_spec(binomial())
  yb <- c(1, 0, 1)
  eb <- c(0.2, -0.1, 0.3)
  qb <- koine:::.re_integrated_loglik(fb, yb, eb, tau = 1.1,
                                                aux = list(), gh = gh)
  bb <- vapply(seq_along(yb), function(i) {
    log(stats::integrate(function(b) {
      p <- stats::plogis(eb[i] + b)
      p^yb[i] * (1 - p)^(1 - yb[i]) * stats::dnorm(b, 0, 1.1)
    }, -Inf, Inf)$value)
  }, numeric(1))
  expect_equal(qb, bb, tolerance = 1e-6)

  fp <- as_family_spec(poisson())
  qp <- koine:::.re_integrated_loglik(fp, 3, 0.5, tau = 0.8,
                                                aux = list(), gh = gh)
  bp <- log(stats::integrate(
    function(b) stats::dpois(3, exp(0.5 + b)) * stats::dnorm(b, 0, 0.8),
    -Inf, Inf)$value)
  expect_equal(qp, bp, tolerance = 1e-4)
})

test_that("Gaussian group marginal matches dense multivariate-normal density", {
  set.seed(5)
  sigma <- 0.7
  tau <- 1.3
  grp <- rep(1:3, c(2, 3, 4))
  y <- stats::rnorm(9)
  eta <- stats::rnorm(9)
  ana <- koine:::.re_group_marginal_gaussian(y, eta, grp, tau, sigma)
  dense <- sum(vapply(1:3, function(g) {
    idx <- grp == g
    m <- sum(idx)
    cov <- sigma^2 * diag(m) + tau^2 * matrix(1, m, m)
    r <- (y - eta)[idx]
    as.numeric(-0.5 * (m * log(2 * pi) +
                         determinant(cov, logarithm = TRUE)$modulus +
                         t(r) %*% solve(cov) %*% r))
  }, numeric(1)))
  expect_equal(ana, dense)
})

test_that("binomial group marginal matches direct integration", {
  gh <- koine:::.gh_rule(24L)
  set.seed(3)
  fb <- as_family_spec(binomial())
  yg <- rbinom(4, 1, 0.6)
  eg <- c(0.1, -0.2, 0.3, 0)
  q <- koine:::.re_group_marginal_loglik(fb, yg, eg, rep(1L, 4), tau = 1.1,
                                         aux = list(), gh = gh)
  brute <- log(stats::integrate(function(b) {
    vapply(b, function(bb) prod(stats::dbinom(yg, 1, stats::plogis(eg + bb))),
           numeric(1)) * stats::dnorm(b, 0, 1.1)
  }, -Inf, Inf)$value)
  expect_equal(q, brute, tolerance = 1e-5)
})

test_that("as tau -> 0 the integrated predictive collapses to the conditional", {
  fb <- as_family_spec(binomial())
  q0 <- koine:::.re_integrated_loglik(fb, 1, 0.2, tau = 1e-8,
                                                aux = list())
  expect_equal(q0, fb$loglik_pointwise(1, 0.2, list()), tolerance = 1e-6)
})
