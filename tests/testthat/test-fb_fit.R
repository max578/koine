make_prov <- function(...) {
  provenance_record(
    method = "dram", method_version = "0.0.0.9000", impl_language = "R",
    seed = 1L, derived_seed = 2L, ir_hash = "abc",
    maturity_tier = "reference-validated",
    independence_key = list(family = "random_walk", kernel = "dram"),
    role = "opinion", exactness = "exact", ...)
}

test_that("provenance validation catches an incomplete tuple", {
  expect_error(
    provenance_record(method = "dram", method_version = "0.1",
                      impl_language = "R", seed = 1L, derived_seed = 2L,
                      ir_hash = "abc", maturity_tier = "reference-validated",
                      independence_key = list(family = "random_walk"),
                      role = "opinion", exactness = "exact"),
    class = "koine_error_provenance_incomplete")
})

test_that("a contributed fb_fit requires draws and validates", {
  arr <- posterior::as_draws_array(
    array(rnorm(2 * 3 * 2), dim = c(3, 2, 2),
          dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
  fit <- fb_fit(draws = arr, provenance = make_prov(), exactness = "exact")
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$status, "contributed")
  expect_null(fit$sample_stats)          # honest NULL (F2)
  expect_null(fit$log_likelihood)        # unassembled (F3)
})

test_that("exactness must agree between fit and provenance", {
  arr <- posterior::as_draws_array(
    array(rnorm(3 * 2 * 1), dim = c(3, 2, 1), dimnames = list(NULL, NULL, "b_x")))
  expect_error(
    fb_fit(draws = arr, provenance = make_prov(), exactness = "resampled"),
    "exactness disagrees")
})

test_that("abstention is first-class and reasons are closed", {
  ab <- abstention("no_diverse_battery", detail = "only one lens")
  expect_true(is_abstention(ab))
  expect_identical(ab$reason, "no_diverse_battery")
  expect_error(abstention("not_a_reason", detail = "x"), "must be one of")
  expect_error(abstention("no_diverse_battery"), "`detail` is required")
})

test_that("an abstained fb_fit carries the abstention, no draws needed", {
  ab <- abstention("loglik_unavailable", scope = "lens", detail = "no ll")
  fit <- fb_fit(draws = NULL, provenance = make_prov(), exactness = "exact",
                status = "abstained", abstention = ab)
  expect_identical(fit$status, "abstained")
  expect_true(is_abstention(fit$abstention))
})
