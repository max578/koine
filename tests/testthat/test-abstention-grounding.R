# FX-11: fact-grounding abstention reasons (Independent Oracle Principle).
# koine can now decline on a FACT failure (un-grounded input), not only a
# numerical one (convergence / calibration / ESS).

test_that("the fact-grounding abstention reasons are in the closed taxonomy", {
  a <- abstention("data_ungrounded",
                  detail = "input NDVI series carries an ungrounded provenance tag")
  expect_true(is_abstention(a))
  expect_identical(a$reason, "data_ungrounded")

  b <- abstention("provenance_unverified", scope = "lens",
                  detail = "upstream manifest obs_schema not grounded")
  expect_true(is_abstention(b))
  expect_identical(b$reason, "provenance_unverified")
})

test_that("an unknown abstention reason is still rejected (closed taxonomy)", {
  expect_error(abstention("made_up_reason", detail = "x"),
               "must be one of")
})

test_that("the numerical abstention reasons still work (no regression)", {
  expect_true(is_abstention(
    abstention("convergence_failed", detail = "rhat > 1.1")))
})
