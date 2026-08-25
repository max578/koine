# K3: DESCRIPTION declares `License: MIT + file LICENSE` and LICENSE/LICENSE.md
# carry the MIT text; README.md and CITATION.cff must not contradict this by
# claiming no licence is declared.

test_that("README does not contradict the declared MIT licence", {
  readme <- readLines(system.file("..", "README.md", package = "koine",
                                  mustWork = FALSE))
  if (length(readme) == 0L) {
    readme <- readLines(file.path(testthat::test_path(), "..", "..",
                                  "README.md"))
  }
  txt <- paste(readme, collapse = "\n")
  expect_false(grepl("No licence is declared yet", txt, fixed = TRUE))
  expect_false(grepl("not yet distributed", txt, fixed = TRUE))
  expect_true(grepl("MIT", txt, fixed = TRUE))
})

test_that("CITATION.cff declares the MIT licence matching DESCRIPTION", {
  cff <- readLines(file.path(testthat::test_path(), "..", "..",
                             "CITATION.cff"))
  txt <- paste(cff, collapse = "\n")
  expect_false(grepl("No licence is declared", txt, fixed = TRUE))
  expect_true(grepl("^license:\\s*MIT\\s*$", txt, perl = TRUE) ||
                any(grepl("^license:\\s*MIT\\s*$", cff, perl = TRUE)))
  desc_license <- read.dcf(file.path(testthat::test_path(), "..", "..",
                                     "DESCRIPTION"))[, "License"]
  expect_true(grepl("MIT", desc_license, fixed = TRUE))
})
