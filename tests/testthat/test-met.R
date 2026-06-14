# Stage-3 MET (multi-environment-trial) structure: two crossed/nested random
# intercepts (variety + trial). Tests the parse, the exact crossed Gaussian
# marginal (validated against an independent lme4 oracle), the integrated
# predictive, the MET fitting path, the registry gates, and back-compatibility
# of the single-random-intercept path.

# A balanced crossed MET simulator with known variance components, used across
# the recovery / corroboration tests (fixed seed -> reproducible).
.sim_met <- function(nv = 12L, nt = 8L, tau_v = 1.3, tau_t = 0.9, sigma = 0.6,
                     beta0 = 5, seed = 101L) {
  set.seed(seed)
  grid <- expand.grid(variety = factor(seq_len(nv)),
                      trial = factor(seq_len(nt)))
  u <- stats::rnorm(nv, 0, tau_v)
  v <- stats::rnorm(nt, 0, tau_t)
  grid$yield <- beta0 + u[grid$variety] + v[grid$trial] +
    stats::rnorm(nrow(grid), 0, sigma)
  attr(grid, "truth") <- list(tau_v = tau_v, tau_t = tau_t, sigma = sigma,
                              beta0 = beta0)
  grid
}

# An UNBALANCED crossed MET (incomplete blocks, each variety in a biased subset
# of trials) with a dominant trial effect --- the realistic case where dropping
# the trial structure confounds the variety effect with which trials each variety
# happened to be grown in. Fixed seed -> reproducible.
.sim_met_unbalanced <- function(nv = 20L, nt = 10L, tau_v = 1.0, tau_t = 3.0,
                                sigma = 0.6, beta0 = 5, seed = 404L) {
  set.seed(seed)
  u <- stats::rnorm(nv, 0, tau_v)
  v <- stats::rnorm(nt, 0, tau_t)
  rows <- list()
  for (i in seq_len(nv)) {
    trials_i <- sample(seq_len(nt), sample(4:7, 1L))
    for (t in trials_i) {
      rows[[length(rows) + 1L]] <- data.frame(
        variety = i, trial = t,
        yield = beta0 + u[i] + v[t] + stats::rnorm(1L, 0, sigma))
    }
  }
  d <- do.call(rbind, rows)
  d$variety <- factor(d$variety)
  d$trial <- factor(d$trial)
  d
}

test_that("two crossed random intercepts parse and flag supported (MET)", {
  d <- .sim_met(nv = 6L, nt = 4L)
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian())
  expect_true(pass$has_re)
  expect_true(pass$re_supported)
  expect_identical(pass$re_structure, "met")
  expect_length(pass$re_terms, 2L)
  red <- .pass_re_design(pass)
  expect_identical(red$structure, "crossed")
  expect_identical(red$groupings, c("variety", "trial"))
  expect_equal(red$factors[[1L]]$n_groups, 6L)
  expect_equal(red$factors[[2L]]$n_groups, 4L)
  # the single-RE legacy fields still mirror the FIRST grouping (back-compat)
  expect_equal(red$n_groups, 6L)
})

test_that("nested groupings are detected as nested, not crossed", {
  # trials nested within regions: each trial sits in exactly one region.
  set.seed(3)
  region <- factor(rep(c("North", "South"), each = 30L))
  trial <- factor(paste0(region, "_", rep(rep(1:3, each = 10L), 2L)))
  d <- data.frame(y = rnorm(60), region = region, trial = trial)
  pass <- koine_pass(y ~ 1 + (1 | region) + (1 | trial), d, gaussian())
  expect_true(pass$re_supported)
  expect_identical(.pass_re_design(pass)$structure, "nested")
})

test_that("more than two groupings, or a random slope, stays unsupported", {
  d <- data.frame(y = rnorm(60), x = rnorm(60),
                  a = factor(rep(1:6, 10)), b = factor(rep(1:5, 12)),
                  c = factor(rep(1:4, 15)))
  expect_false(koine_pass(y ~ x + (1 | a) + (1 | b) + (1 | c), d,
                          gaussian())$re_supported)
  expect_false(koine_pass(y ~ x + (x | a) + (1 | b), d, gaussian())$re_supported)
  # a duplicated grouping is not a valid two-factor MET
  expect_false(koine_pass(y ~ x + (1 | a) + (1 | a), d, gaussian())$re_supported)
})

test_that("the exact crossed marginal equals the lme4 log-likelihood (oracle)", {
  skip_if_not_installed("lme4")
  d <- .sim_met()
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian())
  red <- .pass_re_design(pass)
  design <- .pass_design(pass)
  # independent oracle: lme4 ML fit supplies reference variance components AND a
  # reference marginal log-likelihood koine did not author.
  m <- lme4::lmer(yield ~ 1 + (1 | variety) + (1 | trial), d, REML = FALSE)
  vc <- as.data.frame(lme4::VarCorr(m))
  sigma <- sqrt(vc$vcov[vc$grp == "Residual"])
  tau_v <- sqrt(vc$vcov[vc$grp == "variety"])
  tau_t <- sqrt(vc$vcov[vc$grp == "trial"])
  beta0 <- unname(lme4::fixef(m))
  eta <- rep(beta0, design$N)
  koine_ll <- .re_marginal_gaussian_met(design$y, eta, red$factors,
                                        c(tau_v, tau_t), sigma)
  # the densities are the SAME object -- agreement to machine precision, not a
  # tolerance band (the Independent Oracle Principle satisfied exactly).
  expect_equal(koine_ll, as.numeric(stats::logLik(m)), tolerance = 1e-8)
})

test_that("the crossed marginal reduces to the exact single-RE marginal", {
  # with the second variance component set to zero, the crossed marginal must
  # equal the single-RE per-group Gaussian marginal (a self-consistency anchor).
  d <- .sim_met(nv = 8L, nt = 5L)
  pass2 <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian())
  red2 <- .pass_re_design(pass2)
  design <- .pass_design(pass2)
  eta <- rep(mean(design$y), design$N)
  crossed_one <- .re_marginal_gaussian_met(design$y, eta, red2$factors,
                                           c(1.1, 0), 0.7)
  pass1 <- koine_pass(yield ~ 1 + (1 | variety), d, gaussian())
  red1 <- .pass_re_design(pass1)
  single <- .re_group_marginal_gaussian(design$y, eta, red1$group, 1.1, 0.7)
  expect_equal(crossed_one, single, tolerance = 1e-8)
})

test_that("the MET integrated predictive marginalises every random factor", {
  # a new observation draws a fresh effect from each grouping: the predictive SD
  # is sqrt(sigma^2 + sum tau_r^2).
  y <- c(1, 2, 3)
  eta <- c(0.5, 0.5, 0.5)
  taus <- c(1.2, 0.8)
  sigma <- 0.6
  got <- .re_integrated_loglik_met(y, eta, taus, sigma)
  want <- stats::dnorm(y, eta, sqrt(sigma^2 + sum(taus^2)), log = TRUE)
  expect_equal(got, want)
})

test_that("fit_glmm_marginal recovers a crossed MET posterior", {
  d <- .sim_met()
  truth <- attr(d, "truth")
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian(),
                     seed = 3L)
  fit <- fit_glmm_marginal(pass, n_iter = 2000L, n_burn = 1000L, n_chains = 3L)
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$status, "contributed")
  expect_identical(fit$exactness, "exact")        # the Gaussian MET is exact
  expect_match(fit$provenance$trace$fit, "MET-slice")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "sigma", "tau_variety", "tau_trial"))
  s <- posterior::summarise_draws(fit$draws, "mean", rhat = posterior::rhat)
  m <- stats::setNames(s$mean, s$variable)
  expect_equal(unname(m["sigma"]), truth$sigma, tolerance = 0.2)
  expect_true(m[["tau_variety"]] > 0.4 && m[["tau_variety"]] < 2.5)
  expect_true(m[["tau_trial"]] > 0.3 && m[["tau_trial"]] < 2.5)
  expect_true(all(s$rhat < 1.1))
})

test_that("the MET fit corroborates the production MET model directly", {
  skip_if_not_installed("lme4")
  # an unbalanced MET with a dominant trial effect: the single-RE reduction
  # confounds variety with which trials it was grown in, so it loses fidelity to
  # the production model -- the regime that makes the MET lift worth its cost.
  d <- .sim_met_unbalanced()
  # production MET model (lme4 REML) -> variety BLUPs
  prod <- lme4::lmer(yield ~ 1 + (1 | variety) + (1 | trial), d, REML = TRUE)
  blup_met <- stats::setNames(lme4::ranef(prod)$variety[, 1],
                              rownames(lme4::ranef(prod)$variety))
  # single-RE reduction (variety only, trial structure dropped)
  red1 <- lme4::lmer(yield ~ 1 + (1 | variety), d, REML = TRUE)
  blup_single <- stats::setNames(lme4::ranef(red1)$variety[, 1],
                                 rownames(lme4::ranef(red1)$variety))
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian(),
                     seed = 5L)
  fit <- fit_glmm_marginal(pass, n_iter = 1500L, n_burn = 750L, n_chains = 2L)
  sm <- posterior::summarise_draws(fit$draws, "mean")
  mv <- stats::setNames(sm$mean, sm$variable)
  red <- .pass_re_design(pass)
  design <- .pass_design(pass)
  taus <- c(mv["tau_variety"], mv["tau_trial"])
  v_mat <- .re_marginal_cov(design$N, red$factors, taus, mv["sigma"])
  resid <- design$y -
    as.numeric(design$X %*% mv[paste0("b_", colnames(design$X))])
  z_var <- .re_indicator(red$factors[[1L]]$group, red$factors[[1L]]$n_groups)
  u_koine <- stats::setNames(
    as.numeric(taus[1L]^2 * t(z_var) %*% solve(v_mat, resid)),
    red$factors[[1L]]$levels)
  ord <- names(blup_met)
  r_koine <- stats::cor(u_koine[ord], blup_met[ord])
  r_single <- stats::cor(blup_single[ord], blup_met[ord])
  # koine's MET variety effects corroborate the production MET model tightly,
  # and materially more faithfully than the single-RE reduction (which discards
  # the trial structure the production model and koine both carry).
  expect_gt(r_koine, 0.97)
  expect_gt(r_koine - r_single, 0.1)
})

test_that("the MET marginal concords with the NIMBLE MET reference", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- .sim_met()
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian(),
                     seed = 3L)
  fm <- fit_glmm_marginal(pass, n_iter = 1500L, n_burn = 750L, n_chains = 2L)
  fn <- fit_nimble_glmm(pass, n_iter = 4000L, n_burn = 2000L, n_chains = 2L)
  expect_identical(fn$status, "contributed")
  expect_match(fn$provenance$trace$fit, "MET")
  expect_setequal(posterior::variables(fn$draws),
                  c("b_(Intercept)", "sigma", "tau_variety", "tau_trial"))
  vars <- posterior::variables(fm$draws)
  mm <- stats::setNames(posterior::summarise_draws(fm$draws, "mean")$mean, vars)
  nn <- stats::setNames(posterior::summarise_draws(fn$draws, "mean")$mean,
                        posterior::variables(fn$draws))
  # the marginal (integrates the random effects) and NIMBLE (samples them) carry
  # identical priors, so their close agreement is method-independent consensus.
  expect_lt(abs(mm[["b_(Intercept)"]] - nn[["b_(Intercept)"]]), 0.3)
  expect_lt(abs(mm[["sigma"]] - nn[["sigma"]]), 0.1)
  expect_lt(abs(mm[["tau_variety"]] - nn[["tau_variety"]]), 0.4)
  expect_lt(abs(mm[["tau_trial"]] - nn[["tau_trial"]]), 0.4)
})

test_that("the NIMBLE single-RE lens compiles for an intercept-only fixed part", {
  # regression guard: NIMBLE `inprod` over a length-1 vector generated non-
  # compiling Eigen C++; the intercept-only (K = 1) mean is now scalar. This bit
  # the MET path (its fixed part is commonly intercept-only).
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(4)
  g <- rep(1:10, each = 6)
  d <- data.frame(y = 1 + rnorm(10)[g] + rnorm(60, 0, 0.7), g = factor(g))
  pass <- koine_pass(y ~ 1 + (1 | g), d, gaussian(), seed = 1L)
  fit <- fit_nimble_glmm(pass, n_iter = 1000L, n_burn = 500L, n_chains = 1L)
  expect_identical(fit$status, "contributed")
  expect_setequal(posterior::variables(fit$draws),
                  c("b_(Intercept)", "sigma", "tau"))
})

test_that("a non-Gaussian MET abstains (no exact crossed-GLMM marginal)", {
  set.seed(9)
  grid <- expand.grid(variety = factor(1:10), trial = factor(1:6))
  grid$count <- rpois(nrow(grid), 3)
  pass <- koine_pass(count ~ 1 + (1 | variety) + (1 | trial), grid, poisson(),
                     seed = 1L)
  fit <- fit_glmm_marginal(pass)
  expect_true(is_abstention(fit))
  expect_identical(fit$reason, "feature_unsupported")
  expect_match(fit$detail, "MET")
})

test_that("the MET integrated predictive keeps the PSIS gate reliable", {
  d <- .sim_met(nv = 10L, nt = 6L, seed = 303L)
  pass <- koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d, gaussian(),
                     seed = 4L)
  fit <- fit_glmm_marginal(pass, n_iter = 1200L, n_burn = 600L, n_chains = 2L)
  ll <- assemble_loglik(list(fit), pass)
  expect_equal(ncol(ll[[1L]]$ll), nrow(d))          # per-observation predictive
  gate <- gate_psis(ll)
  expect_length(gate$dropped, 0L)                    # k below the threshold
})

test_that("resolve_battery enrols the MET lenses and excludes single-RE-only ones", {
  d <- .sim_met(nv = 8L, nt = 5L)
  r <- resolve_battery(koine_pass(yield ~ 1 + (1 | variety) + (1 | trial), d,
                                  gaussian()))
  enrolled <- vapply(r$enrolled, function(L) L$id, character(1))
  expect_true("glmm_marginal" %in% enrolled)
  expect_true("nimble_glmm" %in% enrolled)
  # the heavy-tailed single-RE-only lens is excluded with a recorded MET reason
  expect_false(any(grepl("robust_re", enrolled)))
  rre <- Filter(function(e) grepl("robust_re", e$id), r$excluded)
  expect_true(length(rre) >= 1L)
  expect_match(rre[[1L]]$reason, "MET")
  # population lenses stay excluded on a random-effects pass
  expect_false("dram" %in% enrolled)
})

test_that("koine_pass_from_design builds a MET pass from grouping vectors", {
  set.seed(7)
  grid <- expand.grid(variety = factor(1:8), trial = factor(1:5))
  n <- nrow(grid)
  x <- cbind(`(Intercept)` = 1)
  y <- rnorm(n)
  pass <- koine_pass_from_design(
    x[rep(1L, n), , drop = FALSE], y, gaussian(),
    groups = list(variety = grid$variety, trial = grid$trial))
  expect_true(pass$re_supported)
  expect_identical(pass$re_structure, "met")
  red <- .pass_re_design(pass)
  expect_identical(red$structure, "crossed")
  expect_equal(red$factors[[2L]]$n_groups, 5L)
  # supplying both group and groups is an error
  expect_error(
    koine_pass_from_design(x[rep(1L, n), , drop = FALSE], y, gaussian(),
                           group = grid$variety,
                           groups = list(a = grid$variety, b = grid$trial)),
    "OR")
  # more than two groupings is refused
  expect_error(
    koine_pass_from_design(
      x[rep(1L, n), , drop = FALSE], y, gaussian(),
      groups = list(a = grid$variety, b = grid$trial, c = grid$variety)),
    "two random intercepts")
})

test_that("the single random-intercept path is unchanged by the MET lift", {
  # a regression guard: the legacy single-RE accessors and structure tag hold.
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = factor(rep(1:8, 5)))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
  expect_true(pass$re_supported)
  expect_identical(pass$re_structure, "single")
  red <- .pass_re_design(pass)
  expect_identical(red$structure, "single")
  expect_equal(red$n_groups, 8L)
  expect_equal(length(red$group), nrow(d))
})
