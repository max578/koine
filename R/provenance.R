## ----------------------------------------------------------------------------
## Provenance record (Stage-0 SS12, F4)
##
## Mandatory and complete per contributed fb_fit. Cross-engine draw-equality is
## not expected; auditability is. An incomplete provenance tuple on a
## contributed fit is a contract violation (error.provenance_incomplete).
## ----------------------------------------------------------------------------

.provenance_required <- c(
  "method", "method_version", "impl_language", "seed", "derived_seed",
  "ir_hash", "timestamp_iso8601", "maturity_tier", "independence_key",
  "role", "exactness"
)

#' Build a provenance record
#'
#' @param method Lens method id (e.g. `"dram"`).
#' @param method_version Version string of the method implementation.
#' @param impl_language `"R"`, `"C++"`, or `"Python"`.
#' @param seed Master seed from the pass.
#' @param derived_seed Per-lens derived seed.
#' @param ir_hash Content hash of the originating pass.
#' @param maturity_tier One of the maturity tiers (Stage-0 SS11).
#' @param independence_key `list(family =, kernel =)` (Stage-0 SS10).
#' @param role `"opinion"` or `"diagnostic"`.
#' @param exactness `"exact"`, `"approximate"`, or `"resampled"`.
#' @param timestamp_iso8601 ISO-8601 timestamp; defaults to now (UTC).
#' @param lib_versions Optional named list of library versions.
#' @param device Optional device string.
#' @param py_pins Optional named list of Python pins (Tier-3 only).
#' @param trace Optional list of per-stage outcomes.
#' @param structure_variant `TRUE` if the lens fits an alternative model
#'   structure (a different likelihood or prior) rather than the canonical
#'   reference posterior. Carried so the diagnostic layer (the shared-prior
#'   evidence axis, SS9.3) and the independence registry's `specification` axis
#'   (SS10.1) can tell structure variants from reference-posterior lenses.
#' @return A `koine_provenance` object.
#' @examples
#' provenance_record(
#'   method = "dram", method_version = "0.1.0", impl_language = "R",
#'   seed = 1L, derived_seed = 42L, ir_hash = "abc123",
#'   maturity_tier = "reference-validated",
#'   independence_key = list(family = "random_walk", kernel = "dram"),
#'   role = "opinion", exactness = "exact")
#' @export
provenance_record <- function(method, method_version, impl_language, seed,
                              derived_seed, ir_hash, maturity_tier,
                              independence_key, role, exactness,
                              timestamp_iso8601 = NULL, lib_versions = NULL,
                              device = NULL, py_pins = NULL, trace = NULL,
                              structure_variant = FALSE) {
  if (is.null(timestamp_iso8601)) {
    timestamp_iso8601 <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  }
  rec <- list(
    method = method, method_version = method_version,
    impl_language = impl_language, seed = seed, derived_seed = derived_seed,
    ir_hash = ir_hash, timestamp_iso8601 = timestamp_iso8601,
    maturity_tier = maturity_tier, independence_key = independence_key,
    role = role, exactness = exactness, lib_versions = lib_versions,
    device = device, py_pins = py_pins, trace = trace,
    structure_variant = isTRUE(structure_variant)
  )
  class(rec) <- "koine_provenance"
  validate_provenance(rec)
}

#' Validate a provenance record's mandatory tuple
#'
#' @param x A `koine_provenance` object.
#' @return `x`, invisibly, if every mandatory field is present and non-`NULL`;
#'   otherwise raises `error.provenance_incomplete`.
#' @examples
#' prov <- provenance_record(
#'   method = "dram", method_version = "0.1.0", impl_language = "R",
#'   seed = 1L, derived_seed = 42L, ir_hash = "abc123",
#'   maturity_tier = "reference-validated",
#'   independence_key = list(family = "random_walk", kernel = "dram"),
#'   role = "opinion", exactness = "exact")
#' validate_provenance(prov)
#' @export
validate_provenance <- function(x) {
  missing <- .provenance_required[
    vapply(.provenance_required,
           function(f) is.null(x[[f]]) || length(x[[f]]) == 0L,
           logical(1))
  ]
  if (length(missing)) {
    koine_abort("provenance_incomplete",
                paste0("provenance is missing required field(s): ",
                       paste(missing, collapse = ", ")))
  }
  ik <- x$independence_key
  if (!is.list(ik) || is.null(ik$family) || is.null(ik$kernel)) {
    koine_abort("provenance_incomplete",
                "independence_key must be list(family =, kernel =).")
  }
  invisible(x)
}
