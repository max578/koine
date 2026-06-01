test_that("glmm_robust_re recovers a t-RE LMM and emits its own predictive", {
  set.seed(7)
  g <- rep(1:25, each = 8)
  x <- rnorm(200)
  bt <- stats::rt(25, df = 4) * 1.0          # heavy-tailed group effects
  y <- 1 + 0.5 * x + bt[g] + rnorm(200, 0, 0.7)
  d <- data.frame(y = y, x = x, g = factor(g))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian(), seed = 7L)

  fit <- fit_glmm_robust_re(pass, n_iter = 3000L, n_burn = 1500L, n_chains = 3L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$exactness, "approximate")
  expect_false(is.null(fit$log_likelihood))           # native t-RE predictive
  expect_equal(dim(fit$log_likelihood)[3L], nrow(d))
  expect_identical(fit$provenance$independence_key$kernel, "glmm-robust-re-t4")
  expect_true(isTRUE(fit$provenance$structure_variant))

  m <- stats::setNames(posterior::summarise_draws(fit$draws, "mean")$mean,
                       posterior::variables(fit$draws))
  expect_equal(unname(m["b_x"]), 0.5, tolerance = 0.2)
  expect_true(m[["sigma"]] > 0.4 && m[["sigma"]] < 1.0)
})

test_that("glmm_robust_re abstains on a non-Gaussian or population pass", {
  d <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40),
                  g = factor(rep(1:8, 5)))
  expect_true(is_abstention(
    fit_glmm_robust_re(koine_pass(y ~ x + (1 | g), d, binomial()))))
  dp <- data.frame(y = rnorm(30), x = rnorm(30))
  out <- fit_glmm_robust_re(koine_pass(y ~ x, dp, gaussian()))
  expect_true(is_abstention(out))
  expect_identical(out$reason, "feature_unsupported")
})

test_that("the t-RE lens earns more weight on heavy-tailed than Normal effects", {
  # The robust-RE stacking weight when paired against the Gaussian-RE marginal
  # lens. The meaningful, non-flaky claim is RELATIVE --- the t-RE lens earns
  # more weight when the group effects are heavy-tailed than when they are
  # Normal --- so koine's RE opinion is sensitive to random-effect structure.
  robust_weight <- function(seed, re_dist) {
    set.seed(seed)
    g <- rep(1:30, each = 8)
    x <- rnorm(240)
    b <- if (re_dist == "t") {
      stats::rt(30, df = 2.5) * 1.3
    } else {
      rnorm(30, 0, 1.2)
    }
    y <- 1 + 0.5 * x + b[g] + rnorm(240, 0, 0.7)
    pass <- koine_pass(y ~ x + (1 | g), data.frame(y = y, x = x, g = factor(g)),
                       gaussian(), seed = 5L)
    fm <- fit_glmm_marginal(pass, n_iter = 2500L, n_burn = 1200L, n_chains = 2L)
    fr <- fit_glmm_robust_re(pass, n_iter = 2500L, n_burn = 1200L,
                             n_chains = 2L)
    gate <- gate_psis(assemble_loglik(list(fm, fr), pass))
    if (length(gate$kept) < 2L) {
      return(NA_real_)
    }
    w <- stack_opinion(gate$loo)$weights
    unname(w["robust_re_t4"])
  }
  w_t <- robust_weight(22, "t")
  w_normal <- robust_weight(11, "normal")
  skip_if(is.na(w_t) || is.na(w_normal), "a lens dropped at the PSIS gate")
  expect_gt(w_t, w_normal)
})
