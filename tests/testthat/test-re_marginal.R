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

test_that("adaptive GH makes the Poisson group marginal converge across tau", {
  fp <- as_family_spec(poisson())
  yg <- c(1, 7, 1, 2, 9, 4)
  eg <- c(0.1, -0.2, 0.3, 0, 0.5, -0.1)
  for (tau in c(0.5, 1, 1.5, 2, 3, 4)) {
    aghq <- koine:::.re_group_marginal_loglik(fp, yg, eg, rep(1L, 6), tau,
                                              list(), koine:::.gh_rule(20L))
    brute <- log(stats::integrate(function(b) {
      vapply(b, function(z) prod(stats::dpois(yg, exp(eg + z))), numeric(1)) *
        stats::dnorm(b, 0, tau)
    }, -Inf, Inf, rel.tol = 1e-10)$value)
    expect_equal(aghq, brute, tolerance = 1e-3)
  }
})

test_that("marginal gradient matches finite differences (all families)", {
  # the HMC kernel's correctness hinges on this: the analytic marginal gradient
  # (closed form for the Gaussian LMM, the score-identity quadrature estimate
  # for a GLMM) must equal the finite-difference gradient of the value.
  gh <- koine:::.gh_rule(20L)
  fd <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)
  for (fam in c("gaussian", "binomial", "poisson")) {
    set.seed(match(fam, c("gaussian", "binomial", "poisson")))
    fs <- as_family_spec(get(fam, mode = "function")())
    n_g <- 6L
    m <- 5L
    grp <- rep(seq_len(n_g), each = m)
    n <- n_g * m
    eta0 <- stats::rnorm(n, 0, 0.6)
    tau <- 0.8
    sigma <- 1.1
    aux <- if (fam == "gaussian") list(sigma = sigma) else list()
    b <- stats::rnorm(n_g, 0, tau)
    y <- switch(
      fam,
      gaussian = stats::rnorm(n, eta0 + b[grp], sigma),
      binomial = stats::rbinom(n, 1L, stats::plogis(eta0 + b[grp])),
      poisson = stats::rpois(n, exp(eta0 + b[grp])))
    val <- function(eta, tt, ss) {
      a <- if (fam == "gaussian") list(sigma = ss) else list()
      koine:::.re_group_marginal_loglik(fs, y, eta, grp, tt, a, gh)
    }
    g <- koine:::.re_group_marginal_grad(fs, y, eta0, grp, tau, aux, gh)
    # s_obs against d value / d eta_i on a few coordinates
    idx <- c(1L, 12L, 23L, 30L)
    s_fd <- vapply(idx, function(i) {
      fd(function(d) {
        e <- eta0
        e[i] <- e[i] + d
        val(e, tau, sigma)
      }, 0)
    }, numeric(1))
    expect_equal(g$s_obs[idx], s_fd, tolerance = 1e-5)
    # d value / d log tau
    expect_equal(g$d_logtau, fd(function(lt) val(eta0, exp(lt), sigma),
                                log(tau)), tolerance = 1e-5)
    if (fam == "gaussian") {
      expect_equal(g$d_logsigma, fd(function(ls) val(eta0, tau, exp(ls)),
                                    log(sigma)), tolerance = 1e-5)
    } else {
      expect_true(is.na(g$d_logsigma))
    }
  }
})

test_that("Gauss-Laguerre rule integrates against x^alpha exp(-x)", {
  gl <- koine:::.gl_rule(32L, alpha = 1.5)
  expect_equal(sum(gl$weights), gamma(2.5))
  expect_equal(sum(gl$weights * gl$nodes), gamma(3.5))
})

test_that("Student-t scale-mixture marginals match direct integration", {
  nu <- 4
  mix <- koine:::.t_mix_rule(nu, 32L)
  expect_equal(sum(mix$weight), 1)
  expect_equal(sum(mix$weight * mix$lambda), 1)        # Gamma(a, a) mean 1
  tau <- 1.3
  sigma <- 0.7
  # per-observation predictive (b ~ tau * t_nu)
  got <- koine:::.re_integrated_loglik_t(0.4, 0.1, tau, sigma, nu, mix)
  brute <- log(stats::integrate(
    function(b) stats::dnorm(0.4, 0.1 + b, sigma) * stats::dt(b / tau, nu) / tau,
    -Inf, Inf)$value)
  expect_equal(got, brute, tolerance = 1e-4)
  # per-group marginal
  set.seed(1)
  yg <- stats::rnorm(3)
  eg <- c(0.1, -0.2, 0.3)
  gotg <- koine:::.re_group_marginal_t(yg, eg, rep(1L, 3), tau, sigma, nu, mix)
  brug <- log(stats::integrate(function(b) {
    vapply(b, function(z) prod(stats::dnorm(yg, eg + z, sigma)), numeric(1)) *
      stats::dt(b / tau, nu) / tau
  }, -Inf, Inf)$value)
  expect_equal(gotg, brug, tolerance = 1e-5)
})
