# The validation / maturity-gate contract (Stage-0 SS11): SBC, the grid-exact
# reference-posterior check, and the corroboration gate.

test_that("gate_corroboration demotes unbacked koine-first lenses", {
  ref <- gate_corroboration("reference-validated")
  expect_true(ref$may_corroborate)
  expect_equal(ref$discount, 1)

  backed <- gate_corroboration("koine-implemented-validated",
                               validation_ref = "cairn-123")
  expect_true(backed$may_corroborate)
  expect_equal(backed$discount, 0.5)

  unbacked <- gate_corroboration("koine-implemented-validated",
                                 validation_ref = NULL)
  expect_identical(unbacked$effective_tier, "experimental")
  expect_false(unbacked$may_corroborate)
  expect_equal(unbacked$discount, 0)

  exp_tier <- gate_corroboration("experimental")
  expect_false(exp_tier$may_corroborate)
})

test_that("lens_effective_maturity reflects the validation gate", {
  reg <- koine_stage1_registry()
  bid <- function(id) reg[[which(vapply(reg, function(z) z$id, "") == id)]]
  expect_identical(lens_effective_maturity(bid("nimble_gibbs")),
                   "reference-validated")
  # a koine-first lens stripped of its validation_ref is demoted.
  d <- bid("dram")
  d$validation_ref <- NULL
  expect_identical(lens_effective_maturity(d), "experimental")
})

test_that("reference_posterior_check agrees with an analytic 1-D reference", {
  # lens draws ~ N(0, 1); reference log-posterior is the standard normal.
  set.seed(1L)
  arr <- array(stats::rnorm(4000), dim = c(2000L, 2L, 1L),
               dimnames = list(NULL, NULL, "theta"))
  prov <- provenance_record(
    method = "toy", method_version = "0", impl_language = "R",
    seed = 1L, derived_seed = 1L, ir_hash = "h",
    maturity_tier = "reference-validated",
    independence_key = list(family = "x", kernel = "y"),
    role = "opinion", exactness = "exact")
  fit <- fb_fit(draws = posterior::as_draws_array(arr), provenance = prov,
                exactness = "exact", status = "contributed")
  lp <- function(theta) stats::dnorm(theta[["theta"]], 0, 1, log = TRUE)
  chk <- reference_posterior_check(
    fit, lp, grid = list(theta = seq(-6, 6, length.out = 121L)), tol = 0.1)
  expect_true(chk$applicable)
  expect_lt(abs(chk$reference$mean), 0.05)
  expect_lt(abs(chk$reference$sd - 1), 0.05)
  expect_true(chk$passed)
})

test_that("reference_posterior_check refuses an over-large grid", {
  prov <- provenance_record(
    method = "toy", method_version = "0", impl_language = "R",
    seed = 1L, derived_seed = 1L, ir_hash = "h",
    maturity_tier = "reference-validated",
    independence_key = list(family = "x", kernel = "y"),
    role = "opinion", exactness = "exact")
  arr <- array(0, dim = c(2L, 1L, 2L), dimnames = list(NULL, NULL, c("a", "b")))
  fit <- fb_fit(draws = posterior::as_draws_array(arr), provenance = prov,
                exactness = "exact", status = "contributed")
  chk <- reference_posterior_check(
    fit, function(theta) 0,
    grid = list(a = seq_len(100), b = seq_len(100)), max_cells = 100)
  expect_false(chk$applicable)
})

test_that("sbc_check runs the calibration loop and reports uniformity", {
  skip_on_cran()
  setup <- koine_lens_validation_setup(
    "dram", n_obs = 30L, n_pred = 1L,
    fit_control = list(n_iter = 1200L, n_burn = 600L, n_chains = 2L))
  res <- sbc_check(setup$generate, setup$fit_fn, setup$param_names,
                   n_sims = 24L, n_keep = 24L, n_bins = 4L, seed = 7L)
  expect_equal(ncol(res$ranks), length(setup$param_names))
  expect_true(all(res$ranks >= 0 & res$ranks <= res$n_keep))
  expect_true(all(c("param", "chisq", "p_value") %in% names(res$uniformity)))
  expect_type(res$passed, "logical")
})

test_that("reference_posterior_check validates a real DRAM fit on a grid", {
  skip_on_cran()
  setup <- koine_lens_validation_setup(
    "dram", n_obs = 30L, n_pred = 1L,
    fit_control = list(n_iter = 3000L, n_burn = 1500L, n_chains = 3L))
  g <- setup$generate(101L)
  fit <- setup$fit_fn(g$pass)
  ref <- setup$reference(g, n_pts = 25L)
  chk <- reference_posterior_check(fit, ref$log_post_unnorm, ref$grid,
                                   tol = 0.25)
  expect_true(chk$applicable)
  expect_true(all(is.finite(chk$discrepancy)))
})
