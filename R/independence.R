## ----------------------------------------------------------------------------
## The independence-registry contract (Stage-0 SS10)
##
## A two-level key (coarse `family` + fine `kernel`). The grading rule keeps the
## stacking-weight diagnostic from reporting redundant witnesses as cross-axis
## corroboration (the orthogonal-diversity principle).
## ----------------------------------------------------------------------------

# Axis groupings that must collapse to one witness each (Stage-0 SS10). Kernels
# within a grouping are the *same* witness; distinct kernels in the same coarse
# family are discounted, not collapsed.
.koine_axis_collapse <- list(
  population = c("de-mcmc", "dream", "stretch", "t-walk"),
  tempering  = c("static-smc", "drjacoby")
)

#' Build the independence report for the contributing lenses (Stage-0 SS10)
#'
#' One row per contributing opinion lens, carrying its two-level key, maturity
#' tier, and a corroboration discount. The discount reflects: koine-first
#' maturity (SS11), within-family kernel redundancy, and membership of a
#' collapse grouping. This feeds the downstream fourth-vs-three grading.
#'
#' @param fb_fits A list of contributed `koine_fb_fit` objects.
#' @return A data frame with columns `lens_id`, `family`, `kernel`,
#'   `maturity_tier`, `discount`, and `collapsed_witness`.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rnorm(40), x = rnorm(40))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit <- fit_dram(pass, n_iter = 600L, n_burn = 300L, n_chains = 2L)
#' independence_report(list(fit))
#' }
#' @export
independence_report <- function(fb_fits) {
  fits <- Filter(function(f) identical(f$status, "contributed"), fb_fits)
  rows <- lapply(fits, function(f) {
    ik <- f$provenance$independence_key
    grouping <- .collapse_grouping(ik$kernel)
    data.frame(
      lens_id = f$provenance$method,
      family = ik$family,
      kernel = ik$kernel,
      maturity_tier = f$provenance$maturity_tier,
      discount = .maturity_discount(f$provenance$maturity_tier),
      collapsed_witness = if (is.null(grouping)) NA_character_ else grouping,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      lens_id = character(0), family = character(0), kernel = character(0),
      maturity_tier = character(0), discount = numeric(0),
      collapsed_witness = character(0), stringsAsFactors = FALSE)
  }
  out
}

# Maturity discount (Stage-0 SS11): koine-first lenses corroborate at a
# discount; experimental lenses do not corroborate at all.
.maturity_discount <- function(tier) {
  switch(tier,
         "reference-validated" = 1,
         "koine-implemented-validated" = 0.5,
         "experimental" = 0,
         0)
}

# Which collapse grouping (if any) a kernel belongs to.
.collapse_grouping <- function(kernel) {
  for (g in names(.koine_axis_collapse)) {
    if (kernel %in% .koine_axis_collapse[[g]]) return(g)
  }
  NULL
}
