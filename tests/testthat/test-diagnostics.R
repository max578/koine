# The diagnostic layer (Stage-0 SS9): the over-dispersed linear-pool envelope
# and the prior-conditional evidence axis. Both are surfaced, never pooled.

# A minimal contributed fit with given single-chain coefficient draws, so the
# envelope can be exercised without running a sampler.
.fake_fit <- function(method, b_int, b_x, sigma = NULL, kernel = method,
                      structure_variant = FALSE) {
  vars <- list(`b_(Intercept)` = b_int, b_x = b_x)
  if (!is.null(sigma)) vars[["sigma"]] <- sigma
  S <- length(b_int)
  arr <- array(NA_real_, dim = c(S, 1L, length(vars)),
               dimnames = list(NULL, NULL, names(vars)))
  for (j in seq_along(vars)) arr[, 1L, j] <- vars[[j]]
  prov <- provenance_record(
    method = method, method_version = "0", impl_language = "R",
    seed = 1L, derived_seed = 1L, ir_hash = "h",
    maturity_tier = "reference-validated",
    independence_key = list(family = "random_walk", kernel = kernel),
    role = "opinion", exactness = "exact",
    structure_variant = structure_variant)
  fb_fit(draws = posterior::as_draws_array(arr), provenance = prov,
         exactness = "exact", status = "contributed")
}

.toy_gaussian_pass <- function(n = 12L) {
  set.seed(1L)
  d <- data.frame(y = rnorm(n), x = seq_len(n) / n)
  koine_pass(y ~ x, d, gaussian(), seed = 1L)
}

test_that("linear_pool_envelope needs >= 2 comparable lenses, else NULL", {
  pass <- .toy_gaussian_pass()
  one <- list(.fake_fit("a", rnorm(40, 0, 0.1), rnorm(40, 1, 0.1),
                        sigma = rep(1, 40)))
  expect_null(linear_pool_envelope(one, pass))
})

test_that("linear_pool_envelope brackets the pool and flags disagreement", {
  pass <- .toy_gaussian_pass()
  design <- koine:::.pass_design(pass)
  # two lenses agreeing on the slope, one far away on the intercept.
  fits <- list(
    .fake_fit("near_a", rnorm(60, 0, 0.05), rnorm(60, 1, 0.05),
              sigma = rep(1, 60), kernel = "a"),
    .fake_fit("near_b", rnorm(60, 0.02, 0.05), rnorm(60, 1, 0.05),
              sigma = rep(1, 60), kernel = "b"),
    .fake_fit("far", rnorm(60, 5, 0.05), rnorm(60, 1, 0.05),
              sigma = rep(1, 60), kernel = "c"))
  env <- linear_pool_envelope(fits, pass)
  expect_type(env, "list")
  expect_equal(env$n_components, 3L)
  expect_length(env$pool_mean, design$N)
  expect_true(all(env$lower <= env$upper))
  expect_true(all(env$pool_mean >= env$lower & env$pool_mean <= env$upper))
  # the far-off lens disagrees on observations; the near ones largely do not.
  far <- env$flags$n_outside[env$flags$lens_id == "far"]
  near <- env$flags$n_outside[env$flags$lens_id == "near_a"]
  expect_gt(far, near)
})

test_that("a larger inflate widens the envelope", {
  pass <- .toy_gaussian_pass()
  fits <- list(
    .fake_fit("a", rnorm(80, 0, 0.1), rnorm(80, 1, 0.1), sigma = rep(1, 80),
              kernel = "a"),
    .fake_fit("b", rnorm(80, 1, 0.1), rnorm(80, 1, 0.1), sigma = rep(1, 80),
              kernel = "b"))
  narrow <- linear_pool_envelope(fits, pass, inflate = 1)
  wide <- linear_pool_envelope(fits, pass, inflate = 4)
  expect_true(all((wide$upper - wide$lower) >= (narrow$upper - narrow$lower)))
})

test_that("evidence_axis is prior-conditional and excludes structure variants", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")
  pass <- .toy_gaussian_pass(n = 40L)
  ref <- fit_dram(pass, n_iter = 1500L, n_burn = 750L, n_chains = 2L)
  variant <- fit_robust_t(pass, nu = 4, n_iter = 1500L, n_burn = 750L,
                          n_chains = 2L)
  ev <- evidence_axis(list(ref, variant), pass)
  expect_true(ev$prior_conditional)
  expect_identical(ev$role, "diagnostic")
  expect_true(ev$available)
  # the reference-prior lens contributes; the structure variant is excluded.
  expect_true("dram" %in% ev$evidence$lens_id)
  expect_false("robust_t" %in% ev$evidence$lens_id)
  expect_true(is.finite(ev$evidence$log_ml[ev$evidence$lens_id == "dram"]))
})

test_that("the evidence axis never enters log-likelihood assembly (SS6.1)", {
  # a diagnostic-role fit must be skipped by assemble_loglik.
  pass <- .toy_gaussian_pass()
  diag_fit <- .fake_fit("bridge_evidence", rnorm(20), rnorm(20),
                        sigma = rep(1, 20))
  diag_fit$provenance$role <- "diagnostic"
  ll <- assemble_loglik(list(diag_fit), pass)
  expect_length(ll, 0L)
})
