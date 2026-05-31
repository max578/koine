# The diversity engine: a model-structure variant (robust-t) produces a
# genuinely different predictive, so on heavy-tailed data the stacking weights
# do NOT collapse --- koine's orthogonal-diversity premise made visible.

test_that("weight_concentration flags degeneracy by effective component count", {
  even <- weight_concentration(c(a = 0.5, b = 0.5))
  expect_equal(even$n_eff, 2)
  expect_false(even$degenerate)

  conc <- weight_concentration(c(a = 1, b = 0))
  expect_equal(conc$n_eff, 1)
  expect_true(conc$degenerate)

  three <- weight_concentration(c(a = 1 / 3, b = 1 / 3, c = 1 / 3))
  expect_equal(three$n_eff, 3)
})

# small-iter DRAM + robust-t registry
diversity_registry <- function() {
  reg <- koine_stage1_registry()
  bid <- function(id) which(vapply(reg, function(z) z$id, "") == id)
  d <- reg[[bid("dram")]]
  r <- reg[[bid("robust_t")]]
  d$fit_fn <- function(p, ls) fit_dram(p, ls, 2500L, 1200L, 2L)
  r$fit_fn <- function(p, ls) fit_robust_t(p, ls, nu = 4, n_iter = 2500L,
                                           n_burn = 1200L, n_chains = 2L)
  list(d, r)
}

test_that("robust-t realises genuine stacking diversity on heavy-tailed data", {
  skip_on_cran()
  set.seed(42)
  n <- 150L
  x <- rnorm(n)
  y <- 1 + 0.8 * x + stats::rt(n, df = 5) * 0.6      # moderate heavy tails
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 4L)

  res <- synthesise(pass, registry = diversity_registry())
  expect_s3_class(res, "koine_result")
  wc <- res$diagnostics$weight_concentration
  # both lenses carry real weight --- the weights do not collapse to one
  expect_gt(wc$n_eff, 1.3)
  expect_false(wc$degenerate)
  expect_true(all(res$opinion$weights > 0.1))
})

test_that("on_degenerate = 'abstain' yields a weight_degenerate abstention", {
  # two same-posterior DRAM instances corroborate -> degenerate weights
  base <- koine_stage1_registry()[[1L]]
  mk <- function(id, kernel) {
    L <- base
    L$id <- id
    L$independence_key <- list(family = "random_walk", kernel = kernel)
    L$fit_fn <- function(p, ls) fit_dram(p, ls, 1500L, 750L, 2L)
    L
  }
  two_same <- list(mk("dram_a", "dram-a"), mk("dram_b", "dram-b"))
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(60), x = rnorm(60)),
                     gaussian(), seed = 8L)

  out <- synthesise(pass, registry = two_same, on_degenerate = "abstain")
  expect_true(is_abstention(out))
  expect_identical(out$reason, "weight_degenerate")

  # the default ('flag') instead returns the opinion with the degeneracy noted
  flagged <- synthesise(pass, registry = two_same)
  expect_s3_class(flagged, "koine_result")
  expect_true(flagged$diagnostics$weight_concentration$degenerate)
})
