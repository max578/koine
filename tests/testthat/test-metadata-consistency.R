# K3: DESCRIPTION declares `License: MIT + file LICENSE` and LICENSE/LICENSE.md
# carry the MIT text; README.md and CITATION.cff must not contradict this by
# claiming no licence is declared.
#
# README.md and CITATION.cff are top-level source-tree files, not installed
# alongside the package (they carry no `inst/` copy), so under `R CMD check`
# the tests run against an *installed* copy of the package where neither file
# is reachable by a `test_path()`-relative path. These checks are therefore
# source-tree-only and skip cleanly (never error) when the files cannot be
# found relative to the test file -- e.g. under `R CMD check` -- rather than
# failing the whole suite over a check the installed tree cannot perform.

.pkg_root_file <- function(...) {
  file.path(testthat::test_path(), "..", "..", ...)
}

test_that("README does not contradict the declared MIT licence", {
  readme_path <- .pkg_root_file("README.md")
  testthat::skip_if_not(file.exists(readme_path),
                        "README.md not reachable (not a source-tree check).")
  txt <- paste(readLines(readme_path), collapse = "\n")
  expect_false(grepl("No licence is declared yet", txt, fixed = TRUE))
  expect_false(grepl("not yet distributed", txt, fixed = TRUE))
  expect_true(grepl("MIT", txt, fixed = TRUE))
})

test_that("CITATION.cff declares the MIT licence matching DESCRIPTION", {
  cff_path <- .pkg_root_file("CITATION.cff")
  desc_path <- .pkg_root_file("DESCRIPTION")
  testthat::skip_if_not(
    file.exists(cff_path) && file.exists(desc_path),
    "CITATION.cff/DESCRIPTION not reachable (not a source-tree check).")
  cff <- readLines(cff_path)
  txt <- paste(cff, collapse = "\n")
  expect_false(grepl("No licence is declared", txt, fixed = TRUE))
  expect_true(grepl("^license:\\s*MIT\\s*$", txt, perl = TRUE) ||
                any(grepl("^license:\\s*MIT\\s*$", cff, perl = TRUE)))
  desc_license <- read.dcf(desc_path)[, "License"]
  expect_true(grepl("MIT", desc_license, fixed = TRUE))
})
