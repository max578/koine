## ----------------------------------------------------------------------------
## Shrinkage lenses --- sparsity model-structure variants (diversity engine)
##
## A family of Gaussian linear models distinguished only by the PRIOR on the
## coefficients, fitted by bayesreg's Gibbs sampler (Makalic & Schmidt 2016).
## The likelihood is Gaussian throughout --- koine reconstructs each lens's
## predictive via the pass family --- but the priors place different structural
## assumptions on the coefficient vector and so give different posteriors (and
## predictives) whenever the truth has structure:
##
##   ridge      --- a Gaussian (global-scale) prior; mild, near-flat shrinkage.
##   lasso      --- a Laplace prior; moderate, soft-thresholding sparsity.
##   horseshoe  --- a global-local half-Cauchy prior; aggressive sparsity.
##
## Together with the flat-N(0,5) lenses they span a shrinkage spectrum; koine's
## stacking then weighs "how sparse / how shrunk is the signal?" as a
## model-structure question, the second diversity axis after robust-t.
## ----------------------------------------------------------------------------

# Generic bayesreg shrinkage fit, parameterised by the prior. All the family /
# predictor / abstention logic lives here; the public lenses are thin wrappers.
.fit_bayesreg_shrinkage <- function(pass, lens_spec, prior, n_iter, n_burn,
                                    n_chains) {
  if (!requireNamespace("bayesreg", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the bayesreg package is not installed."))
  }
  if (!(pass$family$name == "gaussian" && pass$family$link == "identity")) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("%s serves a continuous (Gaussian) response only.",
                       lens_spec$id)))
  }
  design <- .pass_design(pass)
  if (design$K < 2L) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("%s needs >= 1 predictor beyond the intercept.",
                       lens_spec$id)))
  }
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    set.seed(dseed + c)
    fit <- NULL
    invisible(utils::capture.output(suppressMessages(
      fit <- bayesreg::bayesreg(
        pass$formula, data = pass$data, model = "gaussian", prior = prior,
        n.samples = n_iter, burnin = n_burn, thin = 1L))))
    chains[[c]] <- .bayesreg_chain_matrix(fit, design)
  }
  draws <- .beta_sigma_chains_to_draws(chains, design)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.1.0", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(
      bayesreg = as.character(utils::packageVersion("bayesreg"))),
    trace = list(fit = sprintf("bayesreg (%s) completed", prior)))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

#' Fit the horseshoe (sparsity) lens
#'
#' Gaussian linear model with a horseshoe prior (Carvalho, Polson & Scott 2010)
#' via `bayesreg`. Resolves only for a continuous (Gaussian, identity) response
#' with at least one predictor; otherwise a lens-scope [abstention()].
#'
#' @param pass A `koine_pass` with a Gaussian (identity) family.
#' @param lens_spec The lens's registry record; defaults to the `horseshoe`
#'   record.
#' @param n_iter,n_burn,n_chains Sampler controls.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @references Carvalho, Polson & Scott (2010). The horseshoe estimator for
#'   sparse signals. *Biometrika* 97(2):465-480.
#' @examplesIf requireNamespace("bayesreg", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' d$y <- 1 + 0.8 * d$x1 + rnorm(60)
#' pass <- koine_pass(y ~ x1 + x2, d, gaussian())
#' fit_horseshoe(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_horseshoe <- function(pass, lens_spec = NULL, n_iter = 2000L,
                          n_burn = 1000L, n_chains = 4L) {
  if (is.null(lens_spec)) lens_spec <- .registry_lens("horseshoe")
  .fit_bayesreg_shrinkage(pass, lens_spec, "hs", n_iter, n_burn, n_chains)
}

#' Fit the Bayesian lasso lens
#'
#' Gaussian linear model with a Laplace (Bayesian lasso) prior (Park & Casella
#' 2008) via `bayesreg` --- a moderate, soft-thresholding shrinkage between
#' ridge and horseshoe. Same resolution rules as [fit_horseshoe()].
#'
#' @inheritParams fit_horseshoe
#' @param lens_spec The lens's registry record; defaults to the `lasso` record.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @references Park & Casella (2008). The Bayesian lasso. *JASA*
#'   103(482):681-686.
#' @examplesIf requireNamespace("bayesreg", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' d$y <- 1 + 0.8 * d$x1 + rnorm(60)
#' pass <- koine_pass(y ~ x1 + x2, d, gaussian())
#' fit_lasso(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_lasso <- function(pass, lens_spec = NULL, n_iter = 2000L,
                      n_burn = 1000L, n_chains = 4L) {
  if (is.null(lens_spec)) lens_spec <- .registry_lens("lasso")
  .fit_bayesreg_shrinkage(pass, lens_spec, "lasso", n_iter, n_burn, n_chains)
}

#' Fit the Bayesian ridge lens
#'
#' Gaussian linear model with a Gaussian (ridge) prior of estimated global
#' scale via `bayesreg` --- mild, near-flat shrinkage anchoring the low-sparsity
#' end of the shrinkage spectrum. Same resolution rules as [fit_horseshoe()].
#'
#' @inheritParams fit_horseshoe
#' @param lens_spec The lens's registry record; defaults to the `ridge` record.
#' @return A `koine_fb_fit`, or a `koine_abstention`.
#' @examplesIf requireNamespace("bayesreg", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' d$y <- 1 + 0.8 * d$x1 + rnorm(60)
#' pass <- koine_pass(y ~ x1 + x2, d, gaussian())
#' fit_ridge(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_ridge <- function(pass, lens_spec = NULL, n_iter = 2000L,
                      n_burn = 1000L, n_chains = 4L) {
  if (is.null(lens_spec)) lens_spec <- .registry_lens("ridge")
  .fit_bayesreg_shrinkage(pass, lens_spec, "ridge", n_iter, n_burn, n_chains)
}

# One bayesreg fit -> an (S x [b_(Intercept), b_<preds>, sigma]) matrix, with
# coefficient columns aligned to the pass design's column order.
.bayesreg_chain_matrix <- function(fit, design) {
  bnames <- paste0("b_", colnames(design$X))
  S <- ncol(fit$beta)
  m <- matrix(NA_real_, nrow = S, ncol = length(bnames) + 1L,
              dimnames = list(NULL, c(bnames, "sigma")))
  m[, "b_(Intercept)"] <- as.numeric(fit$beta0)
  for (nm in rownames(fit$beta)) {
    col <- paste0("b_", nm)
    if (col %in% bnames) m[, col] <- as.numeric(fit$beta[nm, ])
  }
  m[, "sigma"] <- sqrt(as.numeric(fit$sigma2))
  m
}

# Stack per-chain (S x V) matrices (with a sigma column) into a draws_array.
.beta_sigma_chains_to_draws <- function(chains, design) {
  vnames <- colnames(chains[[1L]])
  n_keep <- nrow(chains[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(chains), length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (c in seq_along(chains)) arr[, c, ] <- chains[[c]]
  posterior::as_draws_array(arr)
}
