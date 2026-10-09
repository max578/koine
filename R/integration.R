## ----------------------------------------------------------------------------
## Integration surface --- the contract a parent panel (flexyBayes) drives
##
## koine is human-directly-unusable infrastructure: a parent triangulation
## package hands it a model and a required output format, and koine returns a
## synthesised opinion or a structured abstention. Two flexyBayes-AGNOSTIC
## helpers make that handshake clean (koine knows nothing about flexyBayes's
## classes; the parent owns the IR translation + the backend wrapper):
##
##   * `koine_can_contribute()` --- a cheap pre-flight predicate (no fitting):
##     can koine field a diverse opinion for this model, or must it abstain?
##     The parent's capability gate is then a one-liner.
##   * `koine_opinion_draws()` --- materialise koine's opinion as a posterior
##     `draws_array` over the parameters its contributing lenses share, by
##     sampling each lens in proportion to its LOO-stacking weight. This is the
##     representation a parameter-marginal triangulation consumes (a fourth
##     posterior column).
##
## IMPORTANT (charter SS4 caveat). The mixture is over koine's OWN battery of
## lenses --- a legitimate representation of koine's single synthesised opinion.
## It is NOT a pool across the parent's independent engines (that would
## manufacture false confidence and is forbidden). The stacking weights and the
## stacked predictive elpd ride along as attributes so the parent can always
## surface koine's native predictive-space signal alongside the draws.
## ----------------------------------------------------------------------------

#' Can koine contribute a diverse opinion for this model? (pre-flight)
#'
#' A cheap capability predicate --- it resolves the battery but fits nothing ---
#' for a parent panel's backend-availability gate. Returns `TRUE` when at least
#' two independent opinion lenses would enrol, or a single reason string from
#' the abstention taxonomy otherwise (`"feature_unsupported"`,
#' `"no_diverse_battery"`, or `"pass_invalid"` for a malformed model). Because
#' the PSIS reliability gate and weight-degeneracy check can only be evaluated
#' after fitting, a `TRUE` here is necessarily optimistic: koine may still abstain
#' at [synthesise()] time. The predicate is a *floor*, not a guarantee.
#'
#' @param formula,data,family As for [koine_pass()].
#' @param resources A resources list (see [koine_pass()]); defaults to tier 1.
#' @param registry A registry list; defaults to [koine_stage1_registry()].
#' @return `TRUE`, or a single reason-code string.
#' @seealso [koine_opinion_draws()], [synthesise()]
#' @export
#' @examples
#' d <- data.frame(y = rnorm(40), x = rnorm(40))
#' koine_can_contribute(y ~ x, d, gaussian())
koine_can_contribute <- function(formula, data, family,
                                 resources = list(allowed_tiers = 1L),
                                 registry = koine_stage1_registry()) {
  pass <- tryCatch(
    koine_pass(formula, data, family, resources = resources),
    error = function(e) NULL)
  if (is.null(pass)) {
    return("pass_invalid")
  }
  if (isTRUE(pass$has_re) && !isTRUE(pass$re_supported)) {
    return("feature_unsupported")
  }
  resolved <- resolve_battery(pass, registry)
  n_opinion <- sum(vapply(resolved$enrolled,
                          function(L) identical(L$role, "opinion"),
                          logical(1L)))
  if (n_opinion < 2L) {
    return("no_diverse_battery")
  }
  TRUE
}

#' koine's opinion as a posterior draws array (the stacking mixture)
#'
#' Materialises koine's single synthesised opinion as an equal-weight
#' `posterior::draws_array` over the parameters its contributing lenses share, by
#' drawing from each retained lens in proportion to its LOO-stacking weight. This
#' is the parameter-space representation a parent panel's parameter-marginal
#' triangulation consumes. The stacking weights, the stacked predictive elpd, any
#' per-lens parameters dropped from the shared set, and any scale parameters whose
#' lens set was restricted ride along as attributes (`koine_weights`,
#' `koine_elpd_stacked`, `koine_dropped_vars`, `koine_scale_restricted`).
#'
#' The mixture is over koine's own battery (a legitimate single opinion), never a
#' pool across the parent's independent engines (charter SS4). Each parameter is
#' mixed only over the lenses for which it is comparable: the linear-predictor
#' coefficients (`b_<col>`) over every contributing lens, but a scale parameter
#' (`sigma`, `tau`) only over the family-likelihood lenses --- a likelihood-variant
#' lens (the robust-t lens, which carries its own Student-t likelihood) reports its
#' scale on a different metric, so it would bias the synthesised value if averaged
#' in. The returned draws are therefore a per-parameter marginal representation
#' (each column is its comparable-lens stacking mixture); the joint is approximate,
#' which is exactly what a parameter-marginal triangulation needs.
#'
#' @param result A `koine_result` produced with `synthesise(pass,
#'   keep_fits = TRUE)`. A `koine_abstention` is returned unchanged (koine
#'   contributed no opinion to represent).
#' @param n_draws Total draws in the returned mixture. Default 4000.
#' @param seed Optional integer seed for the mixture sampling; defaults to a seed
#'   derived from the pass's content hash, so the representation is reproducible.
#' @return A one-chain `posterior::draws_array`, or the input abstention.
#' @seealso [synthesise()], [koine_can_contribute()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = 1 + rnorm(80), x = rnorm(80))
#' pass <- koine_pass(y ~ x, d, gaussian())
#' res <- synthesise(pass, keep_fits = TRUE)
#' if (!is_abstention(res)) {
#'   dr <- koine_opinion_draws(res)
#'   posterior::summarise_draws(dr, "mean", "sd")
#' }
#' }
koine_opinion_draws <- function(result, n_draws = 4000L, seed = NULL) {
  if (is_abstention(result)) {
    return(result)
  }
  stopifnot(inherits(result, "koine_result"))
  if (is.null(result$fits)) {
    stop("koine_opinion_draws(): no per-lens fits on the result; ",
         "call synthesise(pass, keep_fits = TRUE).", call. = FALSE)
  }
  weights <- result$opinion$weights                # named simplex over kept ids
  by_id <- list()
  family_lik <- list()        # id -> TRUE when the lens fits the pass family
  for (f in result$fits) {
    by_id[[f$provenance$method]] <- f$draws
    # a lens carrying its OWN pointwise log-likelihood predicts under a different
    # likelihood (e.g. the robust-t lens's Student-t), so its scale parameters
    # (sigma, tau) are on a different scale than the family's -- only its
    # linear-predictor coefficients are comparable.
    family_lik[[f$provenance$method]] <- is.null(f$log_likelihood)
  }
  ids <- names(weights)
  if (!all(ids %in% names(by_id))) {
    stop("koine_opinion_draws(): a weighted lens has no retained fit.",
         call. = FALSE)
  }
  mats <- lapply(ids, function(id) {
    unclass(posterior::as_draws_matrix(by_id[[id]]))
  })
  names(mats) <- ids
  shared <- Reduce(intersect, lapply(mats, colnames))
  if (length(shared) == 0L) {
    stop("koine_opinion_draws(): contributing lenses share no common ",
         "parameter to mix over.", call. = FALSE)
  }
  if (is.null(seed)) {
    seed <- as.integer(strtoi(substr(result$ir_hash, 1L, 7L), base = 16L))
  }
  set.seed(seed)
  # Mix each parameter only over the lenses for which it is comparable: the
  # linear-predictor coefficients (`b_<col>`) over every contributing lens; a
  # scale parameter (`sigma`, `tau`, ...) only over the family-likelihood lenses,
  # since a likelihood-variant lens reports it on a different scale (a Student-t
  # scale vs a Gaussian SD). This keeps the synthesised value on one meaning --
  # the structure-variant scale is surfaced via the weights, never averaged into
  # the family scale (charter SS4). The draws are thus a per-parameter marginal
  # representation (the column for each parameter is its comparable-lens stacking
  # mixture); the joint is approximate, which is what a parameter-marginal
  # triangulation consumes.
  restricted <- character(0)
  cols <- lapply(shared, function(v) {
    is_coef <- startsWith(v, "b_")
    sub <- if (is_coef) {
      ids
    } else {
      ids[vapply(ids, function(id) isTRUE(family_lik[[id]]), logical(1L))]
    }
    if (!length(sub) || sum(weights[sub]) == 0) {
      sub <- ids                                 # degenerate: fall back to all
    } else if (!is_coef && length(sub) < length(ids)) {
      restricted <<- c(restricted, v)
    }
    alloc_v <- .allocate_counts(weights[sub], n_draws)
    vals <- unlist(lapply(sub, function(id) {
      n_id <- alloc_v[[id]]
      if (n_id == 0L) {
        return(numeric(0))
      }
      col <- mats[[id]][, v]
      col[sample.int(length(col), n_id, replace = length(col) < n_id)]
    }), use.names = FALSE)
    vals[sample.int(length(vals))]               # shuffle within the column
  })
  mix <- do.call(cbind, cols)
  colnames(mix) <- shared
  arr <- posterior::as_draws_array(
    array(mix, dim = c(nrow(mix), 1L, ncol(mix)),
          dimnames = list(NULL, NULL, colnames(mix))))
  attr(arr, "koine_weights") <- weights
  attr(arr, "koine_elpd_stacked") <- result$opinion$elpd_stacked
  attr(arr, "koine_dropped_vars") <-
    setdiff(Reduce(union, lapply(mats, colnames)), shared)
  attr(arr, "koine_scale_restricted") <- unique(restricted)
  arr
}

# Allocate `n` integer draws across simplex weights `w` by largest remainder, so
# the per-lens counts sum exactly to `n` (the dominant weight absorbs rounding).
.allocate_counts <- function(w, n) {
  w <- w / sum(w)
  raw <- w * n
  base <- floor(raw)
  rem <- as.integer(round(n - sum(base)))
  if (rem > 0L) {
    ord <- order(raw - base, decreasing = TRUE)
    base[ord[seq_len(rem)]] <- base[ord[seq_len(rem)]] + 1
  }
  stats::setNames(as.integer(base), names(w))
}
