## ----------------------------------------------------------------------------
## fb_fit --- the uniform lens-output contract (Stage-0 SS4)
##
## One fb_fit per contributing lens, mapped 1:1 onto an ArviZ DataTree group
## layout so the on-disk NetCDF *is* the contract. Field rules F1-F4 are
## enforced in the validator.
## ----------------------------------------------------------------------------

#' Construct an `fb_fit`
#'
#' @param draws A `posterior::draws_array` of equal-weight posterior draws
#'   (F1). Weighted samplers must resample to equal weight before constructing
#'   the fit and set `exactness = "resampled"`.
#' @param provenance A `koine_provenance` record (mandatory, F4).
#' @param exactness `"exact"`, `"approximate"`, or `"resampled"`.
#' @param sample_stats Optional `posterior::draws_array` of HMC-style sample
#'   statistics; honest `NULL` for non-HMC lenses (F2).
#' @param log_likelihood Optional `array[chain, draw, obs_id]`; usually `NULL`
#'   at construction and filled by [assemble_loglik()] (F3).
#' @param posterior_predictive Optional predictive draws array.
#' @param diagnostics Optional named list of convergence diagnostics; if `NULL`,
#'   `rhat`/`ess` are recovered from `draws` via posterior on demand.
#' @param status `"contributed"` or `"abstained"`.
#' @param abstention An `koine_abstention` iff `status == "abstained"`.
#' @return A validated object of class `koine_fb_fit`.
#' @examples
#' set.seed(1)
#' draws <- posterior::as_draws_array(
#'   array(rnorm(8L), dim = c(2L, 2L, 2L),
#'         dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
#' prov <- provenance_record(
#'   method = "demo", method_version = "0", impl_language = "R",
#'   seed = 1L, derived_seed = 1L, ir_hash = "h",
#'   maturity_tier = "reference-validated",
#'   independence_key = list(family = "random_walk", kernel = "demo"),
#'   role = "opinion", exactness = "exact")
#' fb_fit(draws = draws, provenance = prov, exactness = "exact",
#'        status = "contributed")
#' @export
fb_fit <- function(draws, provenance,
                   exactness = c("exact", "approximate", "resampled"),
                   sample_stats = NULL, log_likelihood = NULL,
                   posterior_predictive = NULL, diagnostics = NULL,
                   status = c("contributed", "abstained"),
                   abstention = NULL) {
  exactness <- match.arg(exactness)
  status <- match.arg(status)
  fit <- structure(
    list(draws = draws, sample_stats = sample_stats,
         log_likelihood = log_likelihood,
         posterior_predictive = posterior_predictive,
         diagnostics = diagnostics, provenance = provenance,
         exactness = exactness, status = status, abstention = abstention),
    class = "koine_fb_fit"
  )
  validate_fb_fit(fit)
}

#' Validate an `fb_fit` against the field rules
#'
#' Enforces: a contributed fit carries draws + a complete provenance tuple (F4);
#' `exactness` agrees between the fit and its provenance; an abstained fit
#' carries an `koine_abstention` and no draws; `sample_stats`/`log_likelihood`
#' are either `NULL` or arrays (no fabricated placeholders, I1).
#'
#' @param x A `koine_fb_fit`.
#' @return `x`, invisibly, if valid; otherwise stops.
#' @examples
#' set.seed(1)
#' draws <- posterior::as_draws_array(
#'   array(rnorm(8L), dim = c(2L, 2L, 2L),
#'         dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
#' prov <- provenance_record(
#'   method = "demo", method_version = "0", impl_language = "R",
#'   seed = 1L, derived_seed = 1L, ir_hash = "h",
#'   maturity_tier = "reference-validated",
#'   independence_key = list(family = "random_walk", kernel = "demo"),
#'   role = "opinion", exactness = "exact")
#' validate_fb_fit(fb_fit(draws, prov, "exact", status = "contributed"))
#' @export
validate_fb_fit <- function(x) {
  if (identical(x$status, "abstained")) {
    if (!is_abstention(x$abstention)) {
      stop("validate_fb_fit(): an abstained fit needs a koine_abstention.",
           call. = FALSE)
    }
    return(invisible(x))
  }
  if (!posterior::is_draws(x$draws)) {
    stop("validate_fb_fit(): `draws` must be a posterior draws object.",
         call. = FALSE)
  }
  validate_provenance(x$provenance)
  if (!identical(x$exactness, x$provenance$exactness)) {
    stop("validate_fb_fit(): exactness disagrees with provenance.",
         call. = FALSE)
  }
  if (!is.null(x$sample_stats) && !posterior::is_draws(x$sample_stats)) {
    stop("validate_fb_fit(): `sample_stats` must be NULL or draws.",
         call. = FALSE)
  }
  if (!is.null(x$log_likelihood) && !is.array(x$log_likelihood)) {
    stop("validate_fb_fit(): `log_likelihood` must be NULL or an array.",
         call. = FALSE)
  }
  invisible(x)
}

#' Recover convergence diagnostics from an `fb_fit`'s draws
#'
#' Computes split-Rhat and bulk/tail ESS from the posterior draws via
#' `posterior::`, engine-agnostically. Never fabricates an HMC-only statistic
#' (I1, F2): for non-HMC lenses `sample_stats` stays `NULL` and only the
#' draws-derived diagnostics are returned.
#'
#' @param x A contributed `koine_fb_fit`.
#' @return A named list with `rhat`, `ess_bulk`, `ess_tail` (per variable).
#' @examples
#' set.seed(1)
#' draws <- posterior::as_draws_array(
#'   array(rnorm(80L), dim = c(20L, 2L, 2L),
#'         dimnames = list(NULL, NULL, c("b_(Intercept)", "b_x"))))
#' prov <- provenance_record(
#'   method = "demo", method_version = "0", impl_language = "R",
#'   seed = 1L, derived_seed = 1L, ir_hash = "h",
#'   maturity_tier = "reference-validated",
#'   independence_key = list(family = "random_walk", kernel = "demo"),
#'   role = "opinion", exactness = "exact")
#' fit <- fb_fit(draws = draws, provenance = prov, exactness = "exact",
#'               status = "contributed")
#' fb_fit_diagnostics(fit)
#' @export
fb_fit_diagnostics <- function(x) {
  if (!inherits(x, "koine_fb_fit") || identical(x$status, "abstained")) {
    stop("fb_fit_diagnostics(): need a contributed koine_fb_fit.",
         call. = FALSE)
  }
  s <- posterior::summarise_draws(
    x$draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  list(
    rhat = stats::setNames(s$rhat, s$variable),
    ess_bulk = stats::setNames(s$ess_bulk, s$variable),
    ess_tail = stats::setNames(s$ess_tail, s$variable)
  )
}

#' @export
print.koine_fb_fit <- function(x, ...) {
  cat("<koine fb_fit>\n")
  cat("  status   : ", x$status, "\n", sep = "")
  if (identical(x$status, "abstained")) {
    cat("  abstained: ", x$abstention$reason, " (", x$abstention$detail, ")\n",
        sep = "")
    return(invisible(x))
  }
  cat("  method   : ", x$provenance$method, "  exactness: ", x$exactness,
      "\n", sep = "")
  cat("  family/kernel: ", x$provenance$independence_key$family, " / ",
      x$provenance$independence_key$kernel, "\n", sep = "")
  cat("  draws    : ", paste(dim(posterior::as_draws_array(x$draws)),
                             collapse = " x "),
      " (chain x draw x var)\n", sep = "")
  cat("  log_lik  : ", if (is.null(x$log_likelihood)) "NULL (unassembled)"
      else paste(dim(x$log_likelihood), collapse = " x "), "\n", sep = "")
  cat("  maturity : ", x$provenance$maturity_tier, "\n", sep = "")
  invisible(x)
}
