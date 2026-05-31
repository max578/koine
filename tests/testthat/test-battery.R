# Full binomial battery: DRAM (random_walk) + NIMBLE (gibbs) + PG-Gibbs (gibbs)
# + DREAM (population) all resolve for a logistic model. They target the same
# posterior, so they corroborate (weights may collapse) --- here we assert the
# breadth of witnesses and the independence grading, not the weight values.

test_that("a binomial model resolves the full four-lens battery", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if_not_installed("BayesLogit")
  skip_if_not_installed("BayesianTools")

  set.seed(11)
  n <- 160L
  x <- rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(0.4 + 1.1 * x))
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), binomial(), seed = 3L)

  # small-iter wrappers so the integration runs in reasonable time
  reg <- koine_stage1_registry()
  by_id <- function(id) which(vapply(reg, function(z) z$id, "") == id)
  reg[[by_id("dram")]]$fit_fn <- function(p, ls) fit_dram(p, ls, 1500L, 750L, 2L)
  reg[[by_id("nimble_gibbs")]]$fit_fn <-
    function(p, ls) fit_nimble(p, ls, 1500L, 750L, 2L)
  reg[[by_id("pg_gibbs")]]$fit_fn <-
    function(p, ls) fit_pg_gibbs(p, ls, 1500L, 750L, 2L)
  reg[[by_id("dream_zs")]]$fit_fn <- function(p, ls) fit_dream(p, ls, 8000L)

  resolved <- resolve_battery(pass, reg)
  ids <- vapply(resolved$enrolled, function(L) L$id, character(1))
  expect_setequal(ids, c("dram", "nimble_gibbs", "pg_gibbs", "dream_zs"))

  res <- synthesise(pass, registry = reg)
  expect_s3_class(res, "koine_result")
  expect_equal(res$opinion$n_components, 4L)
  expect_equal(sum(res$opinion$weights), 1, tolerance = 1e-6)
  # three distinct coarse families contributed
  expect_setequal(unique(res$independence$family),
                  c("random_walk", "gibbs", "population"))
  # the two gibbs kernels are distinct (discounted, not collapsed)
  gibbs_kernels <- res$independence$kernel[res$independence$family == "gibbs"]
  expect_setequal(gibbs_kernels, c("nimble-mwg", "pg-exact"))
})
