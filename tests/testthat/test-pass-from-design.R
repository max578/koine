# Conformance of the IR-targetable builder koine_pass_from_design() against the
# parent backend contract's clause 1 (an array-built, re-entrant, deparse-free
# model builder) plus clauses 2 (draws) and 4 (determinism).

test_that("koine_pass_from_design builds a population pass from arrays", {
  set.seed(1)
  z <- rnorm(40)
  x <- cbind(`(Intercept)` = 1, z = z)        # non-syntactic canonical name
  y <- 1 + 0.5 * z + rnorm(40)
  pass <- koine_pass_from_design(x, y, gaussian(), seed = 4L)
  expect_s3_class(pass, "koine_pass")
  expect_identical(pass$built_from, "design")
  expect_false(pass$has_re)
  d <- koine:::.pass_design(pass)
  # the non-syntactic column name survives verbatim (the formula path mangles it)
  expect_identical(colnames(d$X), c("(Intercept)", "z"))
  expect_identical(d$K, 2L)
})

test_that("the array design matches the formula path on the same data", {
  set.seed(2)
  z <- rnorm(50)
  x <- cbind(`(Intercept)` = 1, z = z)
  y <- 0.5 * z + rnorm(50)
  d_arr <- koine:::.pass_design(koine_pass_from_design(x, y, gaussian()))
  d_form <- koine:::.pass_design(
    koine_pass(y ~ z, data.frame(y = y, z = z), gaussian()))
  cols <- c("(Intercept)", "z")
  expect_equal(unname(d_arr$X[, cols]), unname(d_form$X[, cols]))
  expect_equal(d_arr$y, d_form$y)
})

test_that("koine_pass_from_design carries a single random intercept", {
  set.seed(3)
  g <- rep(1:8, each = 5)
  x <- cbind(`(Intercept)` = 1, z = rnorm(40))
  y <- x[, "z"] + rnorm(8)[g] + rnorm(40)
  pass <- koine_pass_from_design(x, y, gaussian(), group = factor(g))
  expect_true(pass$has_re)
  expect_true(pass$re_supported)
  rd <- koine:::.pass_re_design(pass)
  expect_identical(rd$n_groups, 8L)
  expect_length(rd$group, 40L)
})

test_that("the builder is deterministic (clause 4) and re-entrant (clause 1)", {
  set.seed(4)
  x <- cbind(`(Intercept)` = 1, z = rnorm(30))
  y <- x[, "z"] + rnorm(30)
  a <- koine_pass_from_design(x, y, gaussian(), seed = 9L)
  b <- koine_pass_from_design(x, y, gaussian(), seed = 9L)
  expect_identical(a$ir_hash, b$ir_hash)            # same arrays + seed
  # a second, structurally different pass coexists in the same process
  x2 <- cbind(`(Intercept)` = 1, w = rnorm(20))
  c2 <- koine_pass_from_design(x2, rpois(20, 2), poisson(), seed = 1L)
  expect_false(identical(a$ir_hash, c2$ir_hash))
})

test_that("synthesise + koine_opinion_draws flow through an array-built pass", {
  set.seed(5)
  z <- rnorm(80)
  x <- cbind(`(Intercept)` = 1, z = z)
  y <- 1 + 0.5 * z + rnorm(80)
  res <- synthesise(koine_pass_from_design(x, y, gaussian(), seed = 4L),
                    keep_fits = TRUE)
  skip_if(is_abstention(res))
  dr <- koine_opinion_draws(res)                    # clause 2: a draws_array
  expect_s3_class(dr, "draws_array")
  expect_true(all(c("b_(Intercept)", "b_z", "sigma") %in%
                    posterior::variables(dr)))
  # the array path recovers the same posterior as a formula pass (MC-level)
  rf <- synthesise(koine_pass(y ~ z, data.frame(y = y, z = z), gaussian(),
                              seed = 4L), keep_fits = TRUE)
  skip_if(is_abstention(rf))
  ma <- colMeans(posterior::as_draws_matrix(dr))
  mf <- colMeans(posterior::as_draws_matrix(koine_opinion_draws(rf)))
  expect_equal(unname(ma[c("b_(Intercept)", "b_z")]),
               unname(mf[c("b_(Intercept)", "b_z")]), tolerance = 0.05)
})

test_that("an array-built Poisson random-intercept model contributes (enrolled)", {
  # Poisson random-intercept is now enrolled (the funnel-safe blocked kernel +
  # convergence-gated abstention); on an identifiable model the enrolment reaches
  # through the array builder and synthesises an opinion rather than abstaining
  # (it abstained here only while Poisson was held).
  set.seed(6)
  g <- factor(rep(1:10, 6))
  x <- cbind(`(Intercept)` = 1, w = rnorm(60))
  y <- rpois(60, exp(0.3 + 0.2 * x[, "w"]))
  res <- synthesise(koine_pass_from_design(x, y, poisson(), group = g))
  expect_false(is_abstention(res))
})
