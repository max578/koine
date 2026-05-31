## ----------------------------------------------------------------------------
## pass --- the engine-neutral request (Stage-0 SS2)
##
## A single immutable, declarative object. Contains no engine artefact (no Stan
## code, no brms object) --- else koine's independence leaks. This build covers
## the Stage-1 surface: a population-level GLM formula, the data, the family,
## the canonical observation axis, the required outputs, the seed, and the
## resource envelope.
## ----------------------------------------------------------------------------

#' Construct a koine `pass`
#'
#' @param formula A model formula (population-level terms; random effects and
#'   distributional terms are a later stage).
#' @param data A data frame holding every column referenced in `formula`.
#' @param family A base-R `family` object or a known family name; wrapped into
#'   an evaluable [family_spec()].
#' @param required_outputs Character subset of `c("posterior",
#'   "log_likelihood", "posterior_predictive", "sample_stats")`.
#' @param seed Integer master seed; per-lens seeds are derived from it.
#' @param obs_id Optional character vector of observation ids (defaults to
#'   `obs_0001`-style labels of length `nrow(data)`).
#' @param resources A list with `allowed_tiers` (subset of `1:3`) and optional
#'   budgets; defaults to tier 1 only.
#' @param prior Optional structured prior table (unused by the DRAM lens, which
#'   uses its own weakly-informative default; reserved for later lenses).
#' @return A validated object of class `koine_pass`.
#' @export
#' @examples
#' d <- data.frame(y = rnorm(20), x = rnorm(20))
#' koine_pass(y ~ x, d, gaussian())
koine_pass <- function(formula, data, family,
                       required_outputs = c("posterior", "log_likelihood"),
                       seed = 1L, obs_id = NULL,
                       resources = list(allowed_tiers = 1L),
                       prior = NULL) {
  stopifnot(inherits(formula, "formula"), is.data.frame(data))
  n <- nrow(data)
  if (is.null(obs_id)) {
    obs_id <- sprintf("obs_%04d", seq_len(n))
  }
  if (length(obs_id) != n || anyDuplicated(obs_id)) {
    stop("koine_pass(): `obs_id` must be unique and of length nrow(data) (P3).",
         call. = FALSE)
  }
  fam_spec <- as_family_spec(family)
  vars <- all.vars(formula)
  missing_cols <- setdiff(vars, names(data))
  if (length(missing_cols)) {
    stop("koine_pass(): columns named in formula absent from data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  response <- all.vars(formula)[1L]
  pass <- list(
    format_version = "0.1.0",
    formula = formula,
    response = response,
    data = data,
    family = fam_spec,
    prior = prior,
    obs_dim_names = list(obs_id = obs_id, coords = NULL),
    required_outputs = required_outputs,
    seed = as.integer(seed),
    resources = resources
  )
  pass$ir_hash <- .pass_ir_hash(pass)
  class(pass) <- "koine_pass"
  pass
}

# Content hash over a canonical serialisation of the pass (Stage-0 SS2.4),
# excluding ir_hash itself. The hash is the provenance join key.
.pass_ir_hash <- function(pass) {
  canonical <- list(
    format_version = pass$format_version,
    formula = deparse(pass$formula),
    response = pass$response,
    data = pass$data,
    family = c(pass$family$name, pass$family$link),
    obs_id = pass$obs_dim_names$obs_id,
    required_outputs = sort(pass$required_outputs),
    seed = pass$seed,
    allowed_tiers = pass$resources$allowed_tiers
  )
  digest::digest(canonical, algo = "sha256")
}

# Deterministic per-lens seed derived from the master seed and the lens id.
.derived_seed <- function(master_seed, lens_id) {
  as.integer(
    strtoi(substr(digest::digest(list(master_seed, lens_id), algo = "crc32"),
                  1L, 7L), base = 16L)
  )
}

#' Build the design matrix for a pass (population-level terms)
#'
#' @param pass A `koine_pass`.
#' @return A list with `X` (model matrix), `y` (response vector), and `N`.
#' @keywords internal
.pass_design <- function(pass) {
  mf <- stats::model.frame(pass$formula, data = pass$data)
  X <- stats::model.matrix(pass$formula, data = mf)
  y <- stats::model.response(mf)
  list(X = X, y = as.numeric(y), N = nrow(X), K = ncol(X))
}

#' @export
print.koine_pass <- function(x, ...) {
  cat("<koine pass>\n")
  cat("  formula : ", deparse(x$formula), "\n", sep = "")
  cat("  family  : ", x$family$name, " / ", x$family$link, "\n", sep = "")
  cat("  N       : ", length(x$obs_dim_names$obs_id), "\n", sep = "")
  cat("  outputs : ", paste(x$required_outputs, collapse = ", "), "\n",
      sep = "")
  cat("  ir_hash : ", substr(x$ir_hash, 1L, 12L), "...\n", sep = "")
  invisible(x)
}
