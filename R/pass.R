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
# random-effects terms `(lhs | grouping)`. This build supports random *intercept*
# structure --- a single random intercept `(1 | g)` and the breeder-standard
# multi-environment-trial (MET) structure of two crossed or nested random
# intercepts `(1 | variety) + (1 | trial)`. Richer structure (random slopes,
# correlated effects, more than two groupings) is parsed and flagged unsupported
# (`re_supported = FALSE`) so the battery abstains non-silently (I1), never
# mis-modelling it.
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

  # Support gate. Every bar must be a plain random intercept on a single
  # syntactic grouping name; one such bar is the single-RE case, two on distinct
  # groupings is the MET case (crossed or nested). The (crossed vs nested)
  # distinction is data-driven and resolved later from the design, never from
  # the formula --- the formula is identical either way.
  groupings <- vapply(re_terms, function(t) t$grouping, character(1))
  all_plain_intercepts <- has_re &&
    all(vapply(re_terms, function(t) isTRUE(t$random_intercept), logical(1))) &&
    all(grepl("^[A-Za-z.][A-Za-z0-9._]*$", groupings))
  n_re <- length(re_terms)
  re_supported <- all_plain_intercepts && n_re %in% c(1L, 2L) &&
    !anyDuplicated(groupings)
  re_structure <- if (!re_supported) {
    NA_character_
  } else if (n_re == 1L) {
    "single"
  } else {
    "met"                                          # crossed / nested resolved later
  }
  list(fixed_formula = fixed_formula, re_terms = re_terms,
       has_re = has_re, re_supported = re_supported,
       re_structure = re_structure, re_groupings = groupings,
       response = response)
}

#' Construct a koine `pass`
#'
#' @param formula A model formula. Population-level terms, optionally with
#'   `lme4`-style random *intercepts*: a single `(1 | g)`, or the breeder-standard
#'   multi-environment-trial (MET) structure of two crossed or nested random
#'   intercepts `(1 | variety) + (1 | trial)`. Richer random-effects structure
#'   (random slopes, correlated effects, more than two groupings) is parsed and
#'   flagged unsupported so the battery abstains non-silently rather than fitting
#'   the wrong model. Whether two groupings are crossed or nested is resolved
#'   from the data, never from the formula.
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
    re_structure = ir$re_structure,
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
#'   intercept `(1 | group)`; `NULL` for a population model. For the MET
#'   structure of two crossed/nested random intercepts pass `groups` instead.
#' @param groups Optional named list of two grouping vectors (each length
#'   `nrow(x)`) for two crossed/nested random intercepts --- the breeder-standard
#'   MET structure, e.g. `list(variety = ., trial = .)`. The list names become
#'   the grouping-factor names; mutually exclusive with `group`.
#' @param coef_names Optional column names for `x` (defaults to `colnames(x)`).
#' @param response_name Name recorded for the response (default `"y"`).
#' @param grouping_name Name recorded for the (single) grouping factor when
#'   `group` is supplied (default `"g"`); ignored when `groups` is used.
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
koine_pass_from_design <- function(x, y, family, group = NULL, groups = NULL,
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
  if (!is.null(group) && !is.null(groups)) {
    stop("koine_pass_from_design(): supply `group` (single RE) OR `groups` ",
         "(MET), not both.", call. = FALSE)
  }
  # Normalise the single / MET grouping inputs into one list of factors.
  group_list <- .design_group_list(group, groups, grouping_name, n)
  if (is.null(obs_id)) {
    obs_id <- sprintf("obs_%04d", seq_len(n))
  }
  if (length(obs_id) != n || anyDuplicated(obs_id)) {
    stop("koine_pass_from_design(): `obs_id` must be unique, length nrow(x).",
         call. = FALSE)
  }
  fam_spec <- as_family_spec(family)
  design <- list(X = x, y = y, N = n, K = ncol(x))
  has_re <- length(group_list) > 0L
  re_design <- NULL
  re_terms <- list()
  re_structure <- NA_character_
  if (has_re) {
    re_design <- .build_re_design(group_list)
    re_structure <- if (length(group_list) == 1L) "single" else "met"
    re_terms <- lapply(names(group_list), function(nm) {
      list(grouping = nm, lhs = "1", random_intercept = TRUE)
    })
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
  for (nm in names(group_list)) data[[nm]] <- group_list[[nm]]
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
    re_supported = has_re,            # one or two random intercepts, by construction
    re_structure = re_structure,
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

# Normalise the `group` / `groups` design inputs into a named list of factors
# (empty list for a population pass). Validates lengths and the two-grouping MET
# ceiling: more than two groupings is outside this build's supported structure.
.design_group_list <- function(group, groups, grouping_name, n) {
  if (!is.null(group)) {
    if (length(group) != n) {
      stop("koine_pass_from_design(): length(group) must equal nrow(x).",
           call. = FALSE)
    }
    return(stats::setNames(list(factor(group)), grouping_name))
  }
  if (is.null(groups)) {
    return(list())
  }
  if (!is.list(groups) || is.null(names(groups)) || any(!nzchar(names(groups)))) {
    stop("koine_pass_from_design(): `groups` must be a named list of grouping ",
         "vectors.", call. = FALSE)
  }
  if (length(groups) > 2L) {
    stop("koine_pass_from_design(): at most two random intercepts (the MET ",
         "structure) are supported.", call. = FALSE)
  }
  bad <- vapply(groups, function(g) length(g) != n, logical(1))
  if (any(bad)) {
    stop("koine_pass_from_design(): every vector in `groups` must have length ",
         "nrow(x).", call. = FALSE)
  }
  lapply(groups, factor)
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

#' Random-effects design for a pass (single or MET random intercepts)
#'
#' Returns the per-observation group index/indices and group counts for the
#' grouping factor(s) of the random intercept(s). `NULL` for a population pass.
#' Used by the RE lenses (the marginal likelihood) and the integrated predictive.
#'
#' For the single random-intercept case the returned list carries the legacy
#' fields `grouping`, `group`, `levels`, and `n_groups` directly (so the existing
#' single-RE code path is byte-identical). It always also carries `factors`, a
#' per-grouping list of `list(grouping, group, levels, n_groups)`, the `structure`
#' (`"single"` / `"crossed"` / `"nested"`), and the grouping names --- the surface
#' the MET marginal core reads.
#'
#' @param pass A `koine_pass`.
#' @return A list as described above, or `NULL` when the pass has no random
#'   effect.
#' @keywords internal
.pass_re_design <- function(pass) {
  if (!isTRUE(pass$has_re)) {
    return(NULL)
  }
  if (!is.null(pass$re_design)) {
    return(pass$re_design)
  }
  groupings <- vapply(pass$re_terms, function(t) t$grouping, character(1))
  group_list <- stats::setNames(
    lapply(groupings, function(g) factor(pass$data[[g]])), groupings)
  .build_re_design(group_list)
}

# Build the random-effects design from a named list of grouping factors. One
# factor is the single-RE case; two is the MET case, and whether they are crossed
# or nested is resolved here from the data (nested = each level of the second
# grouping appears within exactly one level of the first, or vice versa). The
# single-RE legacy fields are mirrored at the top level so the existing
# single-intercept code path is untouched.
.build_re_design <- function(group_list) {
  factors <- lapply(names(group_list), function(nm) {
    f <- factor(group_list[[nm]])
    list(grouping = nm, group = as.integer(f), levels = levels(f),
         n_groups = nlevels(f))
  })
  names(factors) <- names(group_list)
  structure_kind <- if (length(factors) == 1L) {
    "single"
  } else if (.re_is_nested(factors[[1L]]$group, factors[[2L]]$group)) {
    "nested"
  } else {
    "crossed"
  }
  out <- list(factors = factors, structure = structure_kind,
              groupings = names(group_list))
  # mirror the single-RE fields at the top level (back-compatible accessors)
  out$grouping <- factors[[1L]]$grouping
  out$group <- factors[[1L]]$group
  out$levels <- factors[[1L]]$levels
  out$n_groups <- factors[[1L]]$n_groups
  out
}

# Two groupings are nested when each level of one appears within a single level
# of the other (the breeder case where trials sit inside regions, say). Crossed
# otherwise (each variety grown across many trials, each trial holding many
# varieties --- the usual MET).
.re_is_nested <- function(g1, g2) {
  nested_within <- function(inner, outer) {
    all(tapply(outer, inner, function(o) length(unique(o)) == 1L))
  }
  nested_within(g2, g1) || nested_within(g1, g2)
}

#' @export
print.koine_pass <- function(x, ...) {
  cat("<koine pass>\n")
  cat("  formula : ", deparse(x$formula), "\n", sep = "")
  if (isTRUE(x$has_re)) {
    terms <- vapply(x$re_terms, function(re) {
      sprintf("(%s | %s)", re$lhs, re$grouping)
    }, character(1))
    tag <- if (!isTRUE(x$re_supported)) {
      "  [unsupported structure]"
    } else if (identical(x$re_structure, "met")) {
      "  [MET]"
    } else {
      ""
    }
    cat("  random  : ", paste(terms, collapse = " + "), tag, "\n", sep = "")
  }
  cat("  family  : ", x$family$name, " / ", x$family$link, "\n", sep = "")
  cat("  N       : ", length(x$obs_dim_names$obs_id), "\n", sep = "")
  cat("  outputs : ", paste(x$required_outputs, collapse = ", "), "\n",
      sep = "")
  cat("  ir_hash : ", substr(x$ir_hash, 1L, 12L), "...\n", sep = "")
  invisible(x)
}
