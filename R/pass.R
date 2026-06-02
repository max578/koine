## ----------------------------------------------------------------------------
## pass --- the engine-neutral request (Stage-0 SS2)
##
## A single immutable, declarative object. Contains no engine artefact (no Stan
## code, no brms object) --- else koine's independence leaks. This build covers
## the Stage-1 surface: a population-level GLM formula, the data, the family,
## the canonical observation axis, the required outputs, the seed, and the
## resource envelope.
## ----------------------------------------------------------------------------

# Decompose a model formula into its fixed-effects part and any lme4-style
# random-effects terms `(lhs | grouping)`. This build supports a single random
# intercept on a single grouping factor; richer structure (random slopes,
# multiple or crossed/nested groupings, correlated effects) is parsed and
# flagged unsupported (`re_supported = FALSE`) so the battery abstains
# non-silently (I1), never mis-modelling it.
.parse_formula_ir <- function(formula) {
  response <- all.vars(formula)[1L]
  rhs <- formula[[length(formula)]]
  bars <- list()
  fixed_terms <- character(0)
  collect <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("+"))) {
      collect(expr[[2L]])
      collect(expr[[3L]])
    } else if (is.call(expr) && identical(expr[[1L]], as.name("(")) &&
                 is.call(expr[[2L]]) &&
                 identical(expr[[2L]][[1L]], as.name("|"))) {
      bars[[length(bars) + 1L]] <<- expr[[2L]]
    } else {
      fixed_terms[[length(fixed_terms) + 1L]] <<- deparse1(expr)
    }
  }
  collect(rhs)
  fixed_terms <- fixed_terms[fixed_terms != "1"]   # implicit intercept retained
  fixed_formula <- if (length(fixed_terms)) {
    stats::reformulate(fixed_terms, response = response)
  } else {
    stats::reformulate("1", response = response)
  }
  re_terms <- lapply(bars, function(b) {
    lhs <- deparse1(b[[2L]])
    list(grouping = deparse1(b[[3L]]), lhs = lhs,
         random_intercept = lhs == "1")
  })
  has_re <- length(re_terms) > 0L
  re_supported <- has_re && length(re_terms) == 1L &&
    isTRUE(re_terms[[1L]]$random_intercept) &&
    grepl("^[A-Za-z.][A-Za-z0-9._]*$", re_terms[[1L]]$grouping)
  list(fixed_formula = fixed_formula, re_terms = re_terms,
       has_re = has_re, re_supported = re_supported, response = response)
}

#' Construct a koine `pass`
#'
#' @param formula A model formula. Population-level terms, optionally with a
#'   single `lme4`-style random intercept `(1 | g)` on one grouping factor.
#'   Richer random-effects structure (random slopes, multiple or crossed/nested
#'   groupings, correlated effects) is parsed and flagged unsupported so the
#'   battery abstains non-silently rather than fitting the wrong model.
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
  ir <- .parse_formula_ir(formula)
  pass <- list(
    format_version = "0.1.0",
    formula = formula,
    fixed_formula = ir$fixed_formula,
    re_terms = ir$re_terms,
    has_re = ir$has_re,
    re_supported = ir$re_supported,
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

#' Construct a koine `pass` from design arrays (the IR-targetable builder)
#'
#' The programmatic, re-entrant, deparse-free entry point a parent triangulation
#' panel drives (flexyBayes backend-contract clause 1): build a pass directly
#' from a named design matrix, a response, a family, and an optional grouping
#' vector --- no model formula, no `match.call()` capture, no session-global
#' state. Two passes can be built and synthesised in the same process without
#' interference. The stored design is used verbatim downstream, so non-syntactic
#' canonical column names (e.g. `(Intercept)`) pass through untouched.
#'
#' @param x A numeric design matrix with named columns (the fixed-effects design;
#'   include an `(Intercept)` column if wanted). Column names become the
#'   coefficient names (`b_<col>` in the draws).
#' @param y The response vector (length `nrow(x)`).
#' @param family A base-R `family` object or a known family name.
#' @param group Optional grouping vector (length `nrow(x)`) for a single random
#'   intercept `(1 | group)`; `NULL` for a population model.
#' @param coef_names Optional column names for `x` (defaults to `colnames(x)`).
#' @param response_name Name recorded for the response (default `"y"`).
#' @param grouping_name Name recorded for the grouping factor (default `"g"`).
#' @inheritParams koine_pass
#' @return A validated object of class `koine_pass`.
#' @seealso [koine_pass()], [koine_opinion_draws()]
#' @export
#' @examples
#' set.seed(1)
#' x <- cbind(`(Intercept)` = 1, z = rnorm(40))
#' y <- 1 + 0.5 * x[, "z"] + rnorm(40)
#' pass <- koine_pass_from_design(x, y, gaussian())
#' pass$has_re             # FALSE -- a population model, built without a formula
koine_pass_from_design <- function(x, y, family, group = NULL,
                                   coef_names = colnames(x),
                                   response_name = "y", grouping_name = "g",
                                   required_outputs = c("posterior",
                                                        "log_likelihood"),
                                   seed = 1L, obs_id = NULL,
                                   resources = list(allowed_tiers = 1L),
                                   prior = NULL) {
  x <- as.matrix(x)
  if (is.null(coef_names) || length(coef_names) != ncol(x)) {
    stop("koine_pass_from_design(): `x` must have one name per column ",
         "(set colnames(x) or pass coef_names).", call. = FALSE)
  }
  colnames(x) <- coef_names
  y <- as.numeric(y)
  n <- nrow(x)
  if (length(y) != n) {
    stop("koine_pass_from_design(): length(y) must equal nrow(x).",
         call. = FALSE)
  }
  if (!is.null(group) && length(group) != n) {
    stop("koine_pass_from_design(): length(group) must equal nrow(x).",
         call. = FALSE)
  }
  if (is.null(obs_id)) {
    obs_id <- sprintf("obs_%04d", seq_len(n))
  }
  if (length(obs_id) != n || anyDuplicated(obs_id)) {
    stop("koine_pass_from_design(): `obs_id` must be unique, length nrow(x).",
         call. = FALSE)
  }
  fam_spec <- as_family_spec(family)
  design <- list(X = x, y = y, N = n, K = ncol(x))
  has_re <- !is.null(group)
  re_design <- NULL
  re_terms <- list()
  if (has_re) {
    gf <- factor(group)
    re_design <- list(grouping = grouping_name, group = as.integer(gf),
                      levels = levels(gf), n_groups = nlevels(gf))
    re_terms <- list(list(grouping = grouping_name, lhs = "1",
                          random_intercept = TRUE))
  }
  # The accessors and most lenses use the stored design directly. A few lenses
  # are formula-driven (the bayesreg shrinkage lens), so synthesise an internal
  # formula + frame from the design: the predictor columns under their own names
  # (an all-ones column is dropped -- bayesreg adds its own intercept), with the
  # response. This stays deparse-free and re-entrant -- the formula is built from
  # the arrays, never captured from a front-end call.
  intercept_cols <- which(apply(x, 2L, function(v) all(v == 1)))
  has_intercept <- length(intercept_cols) > 0L
  pred_idx <- setdiff(seq_len(ncol(x)), intercept_cols)
  pred_names <- colnames(x)[pred_idx]
  data <- as.data.frame(x[, pred_idx, drop = FALSE], check.names = FALSE)
  data[[response_name]] <- y
  if (has_re) data[[grouping_name]] <- factor(group)
  formula <- if (length(pred_names)) {
    stats::reformulate(sprintf("`%s`", pred_names), response = response_name,
                       intercept = has_intercept)
  } else {
    stats::reformulate("1", response = response_name)
  }
  pass <- list(
    format_version = "0.1.0",
    built_from = "design",
    formula = formula,
    fixed_formula = formula,
    re_terms = re_terms,
    has_re = has_re,
    re_supported = has_re,            # a single random intercept, by construction
    response = response_name,
    data = data,
    design = design,
    re_design = re_design,
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

#' Build the design matrix for a pass (fixed-effects terms)
#'
#' Uses the pass's fixed-effects formula (random-effects terms stripped), so the
#' population-level design matrix is identical whether or not the pass carries a
#' random intercept.
#'
#' @param pass A `koine_pass`.
#' @return A list with `X` (model matrix), `y` (response vector), `N`, and `K`.
#' @keywords internal
.pass_design <- function(pass) {
  # a design-built pass (koine_pass_from_design) carries the design verbatim;
  # use it directly so non-syntactic column names (e.g. "(Intercept)") survive.
  if (!is.null(pass$design)) {
    return(pass$design)
  }
  ff <- if (!is.null(pass$fixed_formula)) pass$fixed_formula else pass$formula
  mf <- stats::model.frame(ff, data = pass$data)
  X <- stats::model.matrix(ff, data = mf)
  y <- stats::model.response(mf)
  list(X = X, y = as.numeric(y), N = nrow(X), K = ncol(X))
}

#' Random-effects design for a pass (single random-intercept grouping)
#'
#' Returns the per-observation group index and group count for the grouping
#' factor of the random intercept. `NULL` for a population pass. Used by the
#' RE lenses (the per-group marginal likelihood) and the integrated predictive.
#'
#' @param pass A `koine_pass`.
#' @return A list with `grouping`, `group` (integer index), `levels`, and
#'   `n_groups`, or `NULL` when the pass has no random effect.
#' @keywords internal
.pass_re_design <- function(pass) {
  if (!isTRUE(pass$has_re)) {
    return(NULL)
  }
  if (!is.null(pass$re_design)) {
    return(pass$re_design)
  }
  grouping <- pass$re_terms[[1L]]$grouping
  f <- factor(pass$data[[grouping]])
  list(grouping = grouping, group = as.integer(f), levels = levels(f),
       n_groups = nlevels(f))
}

#' @export
print.koine_pass <- function(x, ...) {
  cat("<koine pass>\n")
  cat("  formula : ", deparse(x$formula), "\n", sep = "")
  if (isTRUE(x$has_re)) {
    re <- x$re_terms[[1L]]
    cat("  random  : (", re$lhs, " | ", re$grouping, ")",
        if (!isTRUE(x$re_supported)) "  [unsupported structure]" else "",
        "\n", sep = "")
  }
  cat("  family  : ", x$family$name, " / ", x$family$link, "\n", sep = "")
  cat("  N       : ", length(x$obs_dim_names$obs_id), "\n", sep = "")
  cat("  outputs : ", paste(x$required_outputs, collapse = ", "), "\n",
      sep = "")
  cat("  ir_hash : ", substr(x$ir_hash, 1L, 12L), "...\n", sep = "")
  invisible(x)
}
