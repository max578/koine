test_that("resolve_battery enrols only the implemented, family-fitting lens", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(20), x = rnorm(20)),
                     gaussian())
  res <- resolve_battery(pass)
  ids <- vapply(res$enrolled, function(L) L$id, character(1))
  # the Gaussian-supported implemented lenses all enrol (subset check, so the
  # test survives the battery growing)
  expect_true(all(c("dram", "nimble_gibbs", "dream_zs", "robust_t",
                    "horseshoe", "lasso", "ridge") %in% ids))
  expect_false("pg_gibbs" %in% ids)      # pg_gibbs is binomial-only

  excl <- vapply(res$excluded, function(e) e$id, character(1))
  expect_true("pg_gibbs" %in% excl)
  pg_reason <- res$excluded[[which(excl == "pg_gibbs")]]$reason
  expect_match(pg_reason, "family")
})

test_that("pg_gibbs would resolve for binomial on family grounds", {
  reg <- koine_stage1_registry()
  pg <- reg[[which(vapply(reg, function(z) z$id, "") == "pg_gibbs")]]
  expect_true(pg$family_support(as_family_spec(binomial())))
  expect_false(pg$family_support(as_family_spec(gaussian())))
})

test_that("tier gating excludes lenses outside allowed_tiers", {
  pass <- koine_pass(y ~ x, data.frame(y = rnorm(20), x = rnorm(20)),
                     gaussian(),
                     resources = list(allowed_tiers = integer(0)))
  res <- resolve_battery(pass)
  expect_length(res$enrolled, 0L)
  expect_true(all(vapply(res$excluded,
                         function(e) grepl("tier", e$reason), logical(1))))
})
