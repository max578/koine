test_that("synthesise abstains when fewer than two opinion lenses resolve", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(30), x = rnorm(30)),
                     gaussian())
  # DRAM (resolves) + PG-Gibbs (excluded: binomial-only) for a Gaussian pass
  # must yield a no-diverse-battery abstention that records the excluded lens.
  reg <- koine_stage1_registry()
  by_id <- function(id) reg[[which(vapply(reg, function(z) z$id, "") == id)]]
  one_plus <- list(by_id("dram"), by_id("pg_gibbs"))
  out <- synthesise(pass, registry = one_plus)
  expect_true(is_abstention(out))
  expect_identical(out$reason, "no_diverse_battery")
  dropped_ids <- vapply(out$lenses_dropped, function(e) e$id, character(1))
  expect_true("pg_gibbs" %in% dropped_ids)
})

# A two-lens registry (two DRAM instances with distinct ids/keys) to exercise
# the full assemble -> gate -> stack path end-to-end.
two_lens_registry <- function() {
  base <- koine_stage1_registry()[[1L]]   # the dram record
  mk <- function(id, kernel) {
    L <- base
    L$id <- id
    L$independence_key <- list(family = "random_walk", kernel = kernel)
    L$fit_fn <- function(pass, lens_spec) {
      fit_dram(pass, lens_spec, n_iter = 1200L, n_burn = 600L, n_chains = 2L)
    }
    L
  }
  list(mk("dram_a", "dram-a"), mk("dram_b", "dram-b"))
}

test_that("synthesise produces a stacked opinion with two lenses", {
  set.seed(5)
  n <- 60L
  x <- rnorm(n)
  y <- 0.5 + 1.2 * x + rnorm(n, sd = 0.8)
  pass <- koine_pass(y ~ x, data.frame(y = y, x = x), gaussian(), seed = 21L)

  res <- synthesise(pass, registry = two_lens_registry())
  expect_s3_class(res, "koine_result")
  expect_identical(res$opinion$n_components, 2L)
  w <- res$opinion$weights
  expect_named(w, c("dram_a", "dram_b"))
  expect_equal(sum(w), 1, tolerance = 1e-6)   # a simplex
  expect_true(all(w >= -1e-8))
  expect_true(is.finite(res$opinion$elpd_stacked))
})

test_that("assemble_loglik builds S x N matrices aligned to obs_id", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(25), x = rnorm(25)),
                     gaussian(), seed = 4L)
  fit <- fit_dram(pass, n_iter = 500L, n_burn = 250L, n_chains = 2L)
  ll <- assemble_loglik(list(fit), pass)
  expect_named(ll, "dram")
  m <- ll[["dram"]]$ll
  expect_equal(ncol(m), 25L)
  expect_identical(colnames(m), pass$obs_dim_names$obs_id)
  # spot-check one cell equals the family log-lik at that draw
  dm <- unclass(posterior::as_draws_matrix(fit$draws))
  X <- model.matrix(y ~ x, pass$data)
  beta1 <- as.numeric(dm[1L, c("b_(Intercept)", "b_x")])
  eta1 <- as.numeric(X %*% beta1)
  expect_equal(
    unname(m[1L, 1L]),
    unname(stats::dnorm(pass$data$y[1L], eta1[1L], dm[1L, "sigma"],
                        log = TRUE)))
})

test_that("gate_psis drops nothing for a well-behaved Gaussian fit", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(50), x = rnorm(50)),
                     gaussian(), seed = 8L)
  fit <- fit_dram(pass, n_iter = 1200L, n_burn = 600L, n_chains = 2L)
  ll <- assemble_loglik(list(fit), pass)
  g <- gate_psis(ll)
  expect_length(g$dropped, 0L)
  expect_true(g$threshold <= 0.7)
})
