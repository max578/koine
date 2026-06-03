## ----------------------------------------------------------------------------
## Abstention --- a first-class, honest outcome (Stage-0 SS8)
##
## koine returns an `abstention` (a value, not a condition) whenever a coherent
## fourth opinion cannot be justified. Honest under-claiming beats a wrong
## fourth opinion. The reason set is closed: adding a reason is a contract
## change.
## ----------------------------------------------------------------------------

# The closed abstention-reason taxonomy (Stage-0 SS8).
.koine_abstention_reasons <- c(
  "no_diverse_battery",   # < 2 independent opinion lenses resolved
  "family_unsupported",   # no lens covers pass$family
  "feature_unsupported",  # the IR uses a part no resolved lens supports
  "loglik_unavailable",   # required pointwise log-likelihood not producible
  "psis_unreliable",      # tier-1 gate failed with no viable fallback
  "calibration_failed",   # tier-2 gate; recalibration insufficient
  "resource_exceeded",    # time / memory / tier budget
  "validation_gate",      # a koine-first lens not yet validated
  "weight_degenerate",    # stacking collapsed to one component (diversity lost)
  "convergence_failed"    # a lens's sampler did not converge (poor ess / rhat)
)

#' Construct a structured abstention
#'
#' An abstention records *why* koine declined to issue a fourth opinion, at what
#' scope, and which lenses (if any) were dropped on the way. It is returned as a
#' normal value and surfaces as "abstained" in the downstream triangulation.
#'
#' @param reason Single string from the closed taxonomy (Stage-0 SS8).
#' @param scope Either `"synthesis"` (the whole opinion) or `"lens"` (a single
#'   lens declined).
#' @param detail Human-readable explanation. Required --- an abstention with no
#'   detail is itself dishonest (I1).
#' @param lenses_dropped Optional list of `list(id =, reason =)` records for
#'   lenses dropped before the abstention.
#' @return An object of class `koine_abstention`.
#' @export
#' @examples
#' abstention("no_diverse_battery", detail = "only 1 opinion lens resolved")
abstention <- function(reason,
                       scope = c("synthesis", "lens"),
                       detail,
                       lenses_dropped = NULL) {
  scope <- match.arg(scope)
  if (!is.character(reason) || length(reason) != 1L ||
        !reason %in% .koine_abstention_reasons) {
    stop("abstention(): `reason` must be one of: ",
         paste(.koine_abstention_reasons, collapse = ", "), call. = FALSE)
  }
  if (missing(detail) || !is.character(detail) || length(detail) != 1L ||
        !nzchar(detail)) {
    stop("abstention(): `detail` is required and must be a non-empty string.",
         call. = FALSE)
  }
  structure(
    list(scope = scope, reason = reason, detail = detail,
         lenses_dropped = lenses_dropped),
    class = "koine_abstention"
  )
}

#' Is this object an abstention?
#'
#' @param x Any object.
#' @return `TRUE` if `x` is a `koine_abstention`.
#' @export
is_abstention <- function(x) inherits(x, "koine_abstention")

#' @export
print.koine_abstention <- function(x, ...) {
  cat("<koine abstention>\n")
  cat("  scope : ", x$scope, "\n", sep = "")
  cat("  reason: ", x$reason, "\n", sep = "")
  cat("  detail: ", x$detail, "\n", sep = "")
  if (length(x$lenses_dropped)) {
    cat("  dropped lenses:\n")
    for (d in x$lenses_dropped) {
      cat("    - ", d$id, ": ", d$reason, "\n", sep = "")
    }
  }
  invisible(x)
}
