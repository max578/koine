test_that("a population formula carries no random effect", {
  d <- data.frame(y = rnorm(20), x = rnorm(20))
  pass <- koine_pass(y ~ x, d, gaussian())
  expect_false(isTRUE(pass$has_re))
  expect_equal(.pass_design(pass)$K, 2L)        # intercept + x
  expect_null(.pass_re_design(pass))
})

test_that("a single random intercept is parsed and flagged supported", {
  d <- data.frame(y = rnorm(20), x = rnorm(20), g = factor(rep(1:4, 5)))
  pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
  expect_true(pass$has_re)
  expect_true(pass$re_supported)
  expect_identical(pass$re_terms[[1L]]$grouping, "g")
  expect_equal(.pass_design(pass)$K, 2L)        # fixed design strips the RE
  red <- .pass_re_design(pass)
  expect_equal(red$n_groups, 4L)
  expect_equal(length(red$group), nrow(d))
})

test_that("unsupported random-effects structures are flagged, never mis-parsed", {
  d <- data.frame(y = rnorm(20), x = rnorm(20),
                  g = factor(rep(1:4, 5)), h = factor(rep(1:2, 10)))
  expect_false(koine_pass(y ~ x + (x | g), d, gaussian())$re_supported)  # slope
  expect_false(koine_pass(y ~ x + (1 | g) + (1 | h), d, gaussian())$re_supported)
})

test_that("the feature gate matches lenses to passes symmetrically", {
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = factor(rep(1:8, 5)))
  # population pass: only population (non-RE) lenses enrol
  r0 <- resolve_battery(koine_pass(y ~ x, d, gaussian()))
  enrolled0 <- vapply(r0$enrolled, function(L) L$id, character(1))
  expect_false(any(c("glmm_marginal", "nimble_glmm") %in% enrolled0))
  # RE pass: the population lenses are excluded with a feature reason
  r1 <- resolve_battery(koine_pass(y ~ x + (1 | g), d, gaussian()))
  enrolled1 <- vapply(r1$enrolled, function(L) L$id, character(1))
  expect_false("dram" %in% enrolled1)
  dram_excl <- Filter(function(e) e$id == "dram", r1$excluded)[[1L]]
  expect_match(dram_excl$reason, "random-effects")
})

test_that("synthesise abstains feature_unsupported on an unsupported RE structure", {
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = factor(rep(1:8, 5)))
  res <- synthesise(koine_pass(y ~ x + (x | g), d, gaussian()))
  expect_true(is_abstention(res))
  expect_identical(res$reason, "feature_unsupported")
})
