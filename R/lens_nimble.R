## ----------------------------------------------------------------------------
## NIMBLE-Gibbs lens (Stage-0 SS13; independence family `gibbs / nimble-mwg`)
##
## Wraps NIMBLE's conjugate / Metropolis-within-Gibbs MCMC. It targets the
## *same* posterior as the DRAM lens (identical weakly-informative priors:
## beta ~ N(0, 5); for Gaussian, log_sigma ~ N(0, 2)), so any agreement between
## the two is method-independent consensus and any disagreement is a genuine
## inference-family signal --- the orthogonal-diversity axis (SS2).
##
## NIMBLE compiles generated C++; the lens silences that chatter and runs under
## a deterministic per-lens seed. It returns posterior draws only; koine
## reconstructs the pointwise log-likelihood downstream via the family spec
## (F3), keeping the assembly path uniform across the battery.
## ----------------------------------------------------------------------------

# Build the family-specific NIMBLE model code, or NULL if the family/link pair
# is not supported by this lens (the caller then abstains, lens scope).
.nimble_model_code <- function(family_name, link) {
  ok <- (family_name == "gaussian" && link == "identity") ||
    (family_name == "binomial" && link == "logit") ||
    (family_name == "poisson" && link == "log")
  if (!ok) return(NULL)
  switch(
    family_name,
    gaussian = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = 5)
      log_sigma ~ dnorm(0, sd = 2)
      sigma <- exp(log_sigma)
      for (i in 1:N) {
        mu[i] <- inprod(X[i, 1:K], beta[1:K])
        y[i] ~ dnorm(mu[i], sd = sigma)
      }
    }),
    binomial = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = 5)
      for (i in 1:N) {
        logit(p[i]) <- inprod(X[i, 1:K], beta[1:K])
        y[i] ~ dbin(p[i], 1)
      }
    }),
    poisson = nimble::nimbleCode({
      for (k in 1:K) beta[k] ~ dnorm(0, sd = 5)
      for (i in 1:N) {
        log(lambda[i]) <- inprod(X[i, 1:K], beta[1:K])
        y[i] ~ dpois(lambda[i])
      }
    })
  )
}

# Evaluate an expression with NIMBLE's printed output and messages silenced.
.nimble_quiet <- function(expr) {
  invisible(utils::capture.output(
    suppressMessages(suppressWarnings(force(expr)))))
}

#' Fit the NIMBLE-Gibbs lens
#'
#' Runs NIMBLE's conjugate / Metropolis-within-Gibbs sampler on the pass's GLM
#' and returns an [fb_fit()] of equal-weight posterior draws with a complete
#' provenance record. If the family/link pair is unsupported, returns a
#' lens-scope [abstention()] rather than a wrong fit (I1).
#'
#' @param pass A `koine_pass`.
#' @param lens_spec The lens's registry record; defaults to the registry's
#'   `nimble_gibbs` record.
#' @param n_iter Iterations per chain (including burn-in).
#' @param n_burn Burn-in iterations discarded per chain.
#' @param n_chains Number of chains.
#' @return A `koine_fb_fit`, or a `koine_abstention` for an unsupported link.
#' @examplesIf requireNamespace("nimble", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(40))
#' d$y <- 1 + 0.8 * d$x + rnorm(40)
#' pass <- koine_pass(y ~ x, d, gaussian())
#' fit_nimble(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
#' }
#' @export
fit_nimble <- function(pass, lens_spec = NULL, n_iter = 4000L, n_burn = 2000L,
                       n_chains = 4L) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    return(abstention("resource_exceeded", scope = "lens",
                      detail = "the nimble package is not installed."))
  }
  # NIMBLE's model-build / generated code calls some of its functions
  # unqualified (e.g. getNimbleOption), so nimble must be ON THE SEARCH PATH,
  # not merely namespace-loaded via ::. Attach it on first use.
  if (!"package:nimble" %in% search()) {
    suppressMessages(attachNamespace("nimble"))
  }
  if (is.null(lens_spec)) {
    reg <- koine_stage1_registry()
    lens_spec <- reg[[which(vapply(reg, function(z) z$id, "") ==
                              "nimble_gibbs")]]
  }
  code <- .nimble_model_code(pass$family$name, pass$family$link)
  if (is.null(code)) {
    return(abstention(
      "feature_unsupported", scope = "lens",
      detail = sprintf("nimble_gibbs has no model for family '%s' + link '%s'.",
                       pass$family$name, pass$family$link)))
  }

  design <- .pass_design(pass)
  has_sigma <- identical(pass$family$name, "gaussian")
  dseed <- .derived_seed(pass$seed, lens_spec$id)
  monitors <- if (has_sigma) c("beta", "sigma") else "beta"

  constants <- list(N = design$N, K = design$K, X = design$X)
  data_list <- list(y = design$y)
  warm <- tryCatch(
    suppressWarnings(stats::coef(stats::glm.fit(
      design$X, design$y, family = .pass_glm_family(pass)))),
    error = function(e) rep.int(0, design$K))
  warm[!is.finite(warm)] <- 0
  inits <- list(beta = as.numeric(warm))
  if (has_sigma) inits$log_sigma <- 0

  samples <- NULL
  .nimble_quiet({
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = data_list, inits = inits,
                                 calculate = FALSE)
    cmodel <- nimble::compileNimble(model)
    conf <- nimble::configureMCMC(model, monitors = monitors, print = FALSE)
    mcmc <- nimble::buildMCMC(conf)
    cmcmc <- nimble::compileNimble(mcmc, project = model)
    samples <- nimble::runMCMC(
      cmcmc, niter = n_iter, nburnin = n_burn, nchains = n_chains,
      setSeed = dseed + seq_len(n_chains), progressBar = FALSE,
      samplesAsCodaMCMC = FALSE)
  })
  if (!is.list(samples)) samples <- list(samples)

  draws <- .nimble_to_draws(samples, design, has_sigma)
  prov <- provenance_record(
    method = lens_spec$id, method_version = "0.1.0", impl_language = "R",
    seed = pass$seed, derived_seed = dseed, ir_hash = pass$ir_hash,
    maturity_tier = lens_spec$maturity_tier,
    independence_key = lens_spec$independence_key,
    role = lens_spec$role, exactness = lens_spec$exactness,
    structure_variant = isTRUE(lens_spec$structure_variant),
    lib_versions = list(nimble = as.character(utils::packageVersion("nimble"))),
    trace = list(fit = "nimble MCMC completed"))
  fb_fit(draws = draws, provenance = prov, exactness = lens_spec$exactness,
         sample_stats = NULL, status = "contributed")
}

# Convert NIMBLE's per-chain sample matrices (columns beta[1].. [, sigma]) into
# a posterior draws_array, renaming coefficients to the b_<col> convention
# assemble_loglik expects.
.nimble_to_draws <- function(samples, design, has_sigma) {
  bnames <- paste0("b_", colnames(design$X))
  target <- if (has_sigma) c(bnames, "sigma") else bnames
  n_keep <- nrow(samples[[1L]])
  arr <- array(NA_real_, dim = c(n_keep, length(samples), length(target)),
               dimnames = list(NULL, NULL, target))
  for (c in seq_along(samples)) {
    m <- samples[[c]]
    src <- matrix(NA_real_, nrow = n_keep, ncol = length(target),
                  dimnames = list(NULL, target))
    for (k in seq_len(design$K)) {
      src[, bnames[k]] <- m[, paste0("beta[", k, "]")]
    }
    if (has_sigma) src[, "sigma"] <- m[, "sigma"]
    arr[, c, ] <- src
  }
  posterior::as_draws_array(arr)
}
