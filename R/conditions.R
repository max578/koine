## ----------------------------------------------------------------------------
## Typed conditions: errors vs abstentions (Stage-0 SS8, SS12)
##
## Two disjoint failure channels, never conflated:
##   * an *error* signals a contract violation (a bug, a malformed input, a
##     broken NetCDF file). It halts with a typed condition and is never
##     swallowed (charter hard-rule 2).
##   * an *abstention* is a normal outcome (no diverse battery, an
##     unsupported family, an unreliable PSIS gate). It is a value, not a
##     condition --- see `abstention()`.
## ----------------------------------------------------------------------------

# The closed error taxonomy (Stage-0 SS12). Adding a code is a contract change.
.koine_error_codes <- c(
  "pass_hash_mismatch",
  "provenance_incomplete",
  "netcdf_shape",
  "netcdf_dimorder",
  "obs_axis_mismatch",
  "netcdf_dtype",
  "family_not_evaluable",
  "registry_invalid"
)

#' Signal a koine contract-violation error
#'
#' Raises a classed condition of subclass `koine_error_<code>` (and the parent
#' classes `koine_error` and `error`). The `code` must be one of the closed
#' error taxonomy in Stage-0 SS12; an unknown code is itself an error, because a
#' silent, untyped failure would violate charter hard-rule 2.
#'
#' @param code Single string; one of `.koine_error_codes`.
#' @param message Human-readable message.
#' @param ... Extra named fields attached to the condition for inspection.
#' @return Never returns; always stops.
#' @keywords internal
koine_abort <- function(code, message, ...) {
  if (!is.character(code) || length(code) != 1L ||
        !code %in% .koine_error_codes) {
    stop("internal: koine_abort() called with an unknown error code: ",
         deparse(code), call. = FALSE)
  }
  cond <- structure(
    class = c(paste0("koine_error_", code), "koine_error", "error",
              "condition"),
    list(message = message, call = sys.call(-1L), code = code, ...)
  )
  stop(cond)
}
