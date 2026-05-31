## ----------------------------------------------------------------------------
## family_spec --- the evaluable-likelihood interface (Stage-0 SS2.3, C2)
##
## This is the field that makes synthesis work across a heterogeneous battery.
## Many lenses (DRAM, static SMC, nested sampling, any koine-first sampler)
## yield posterior *draws* but no pointwise log-likelihood matrix. koine
## reconstructs the S x N matrix by evaluating `loglik_pointwise(y, eta, aux)`
## per draw. The upstream family object (a base-R / brms `family`) exposes
## `linkinv` but no pointwise log-likelihood, so koine owns it here.
## ----------------------------------------------------------------------------

#' Construct an evaluable family specification
#'
#' A `family_spec` exposes, as pure functions, the pieces koine needs to compute
#' a pointwise log-likelihood and a posterior-predictive draw for any lens that
#' does not emit them natively. This low-level constructor is normally reached
#' through [as_family_spec()], which adapts a base-R `family` object.
#'
#' @param name Family name, e.g. `"gaussian"`, `"binomial"`, `"poisson"`.
#' @param link Link name, e.g. `"identity"`, `"logit"`, `"log"`.
#' @param linkinv Function mapping the linear predictor `eta` to the
#'   response-scale mean `mu`.
#' @param loglik_pointwise Function `(y, eta, aux)` returning a length-`N`
#'   vector of `log p(y_i | eta_i, aux)`. The contract koine relies on.
#' @param rng Function `(eta, aux)` returning one posterior-predictive draw
#'   `y_tilde` (used only when posterior-predictive output is requested).
#' @param support_check Function `(y)` returning `TRUE` if `y` lies in the
#'   family support.
#' @param aux_params Character vector naming the auxiliary parameters (e.g.
#'   `"sigma"` for Gaussian, none for Binomial).
#' @param dlp_deta Optional function `(y, eta, aux)` returning the first
#'   derivative of the pointwise log-likelihood with respect to the linear
#'   predictor (the canonical-link score), used by adaptive Gauss-Hermite
#'   quadrature to locate the random-effect integrand's mode.
#' @param d2lp_deta2 Optional function `(y, eta, aux)` returning the second
#'   derivative (the negative observed information), used for the quadrature's
#'   curvature scaling.
#' @return An object of class `koine_family_spec`.
#' @seealso [as_family_spec()]
#' @examples
#' # normally reached through as_family_spec(); the low-level constructor:
#' fs <- family_spec(
#'   name = "gaussian", link = "identity",
#'   linkinv = function(eta) eta,
#'   loglik_pointwise = function(y, eta, aux) {
#'     stats::dnorm(y, eta, aux$sigma, log = TRUE)
#'   },
#'   rng = function(eta, aux) stats::rnorm(length(eta), eta, aux$sigma),
#'   support_check = function(y) all(is.finite(y)),
#'   aux_params = "sigma")
#' fs$loglik_pointwise(c(0, 1), c(0, 1), list(sigma = 1))
#' @export
family_spec <- function(name, link, linkinv, loglik_pointwise, rng,
                        support_check, aux_params = character(0),
                        dlp_deta = NULL, d2lp_deta2 = NULL) {
  stopifnot(
    is.character(name), length(name) == 1L,
    is.character(link), length(link) == 1L,
    is.function(linkinv), is.function(loglik_pointwise),
    is.function(rng), is.function(support_check),
    is.character(aux_params)
  )
  structure(
    list(name = name, link = link, aux_params = aux_params,
         linkinv = linkinv, loglik_pointwise = loglik_pointwise,
         rng = rng, support_check = support_check,
         dlp_deta = dlp_deta, d2lp_deta2 = d2lp_deta2),
    class = "koine_family_spec"
  )
}

#' Adapt a base-R family to an evaluable `family_spec`
#'
#' Wraps a base-R / brms `family` object (which exposes `linkinv` but no
#' pointwise log-likelihood) with the `loglik_pointwise` and `rng` functions
#' koine needs. Supported families in this build are `gaussian`, `binomial`
#' (and `bernoulli`, a binomial with one trial), and `poisson`. An unsupported
#' family raises `error.family_not_evaluable`, never a silent approximation.
#'
#' @param family A base-R `family` object (e.g. `gaussian()`, `binomial()`), or
#'   a single string naming one.
#' @param trials Optional length-`N` integer vector of binomial trial counts;
#'   defaults to 1 (the Bernoulli case).
#' @return A `koine_family_spec`.
#' @export
#' @examples
#' fs <- as_family_spec(gaussian())
#' fs$loglik_pointwise(y = c(0, 1), eta = c(0, 1), aux = list(sigma = 1))
as_family_spec <- function(family, trials = NULL) {
  if (is.character(family)) {
    family <- switch(
      family,
      gaussian = stats::gaussian(),
      binomial = stats::binomial(),
      bernoulli = stats::binomial(),
      poisson = stats::poisson(),
      koine_abort("family_not_evaluable",
                  sprintf("unknown family name: '%s'", family))
    )
  }
  if (!inherits(family, "family")) {
    koine_abort("family_not_evaluable",
                "`family` must be a base-R family object or a known name.")
  }
  linkinv <- family$linkinv
  fam <- family$family

  spec <- switch(
    fam,
    gaussian = .family_spec_gaussian(family$link, linkinv),
    binomial = .family_spec_binomial(family$link, linkinv, trials),
    poisson  = .family_spec_poisson(family$link, linkinv),
    koine_abort(
      "family_not_evaluable",
      sprintf("family '%s' has no evaluable spec in this build.", fam))
  )
  spec
}

# -- gaussian ----------------------------------------------------------------
.family_spec_gaussian <- function(link, linkinv) {
  family_spec(
    name = "gaussian", link = link, aux_params = "sigma",
    linkinv = linkinv,
    loglik_pointwise = function(y, eta, aux) {
      stats::dnorm(y, mean = linkinv(eta), sd = aux$sigma, log = TRUE)
    },
    rng = function(eta, aux) {
      stats::rnorm(length(eta), mean = linkinv(eta), sd = aux$sigma)
    },
    support_check = function(y) is.numeric(y) && all(is.finite(y)),
    dlp_deta = function(y, eta, aux) (y - linkinv(eta)) / aux$sigma^2,
    d2lp_deta2 = function(y, eta, aux) rep(-1 / aux$sigma^2, length(y))
  )
}

# -- binomial / bernoulli ----------------------------------------------------
.family_spec_binomial <- function(link, linkinv, trials) {
  family_spec(
    name = "binomial", link = link, aux_params = character(0),
    linkinv = linkinv,
    loglik_pointwise = function(y, eta, aux) {
      n <- if (is.null(trials)) rep.int(1L, length(y)) else trials
      stats::dbinom(y, size = n, prob = linkinv(eta), log = TRUE)
    },
    rng = function(eta, aux) {
      n <- if (is.null(trials)) 1L else trials
      stats::rbinom(length(eta), size = n, prob = linkinv(eta))
    },
    support_check = function(y) {
      all(is.finite(y)) && all(y >= 0) && all(y == round(y))
    },
    # canonical logit-link score / information (n = trials)
    dlp_deta = function(y, eta, aux) {
      n <- if (is.null(trials)) rep.int(1L, length(y)) else trials
      y - n * linkinv(eta)
    },
    d2lp_deta2 = function(y, eta, aux) {
      n <- if (is.null(trials)) rep.int(1L, length(y)) else trials
      p <- linkinv(eta)
      -n * p * (1 - p)
    }
  )
}

# -- poisson -----------------------------------------------------------------
.family_spec_poisson <- function(link, linkinv) {
  family_spec(
    name = "poisson", link = link, aux_params = character(0),
    linkinv = linkinv,
    loglik_pointwise = function(y, eta, aux) {
      stats::dpois(y, lambda = linkinv(eta), log = TRUE)
    },
    rng = function(eta, aux) {
      stats::rpois(length(eta), lambda = linkinv(eta))
    },
    support_check = function(y) {
      all(is.finite(y)) && all(y >= 0) && all(y == round(y))
    },
    # canonical log-link score / information
    dlp_deta = function(y, eta, aux) y - linkinv(eta),
    d2lp_deta2 = function(y, eta, aux) -linkinv(eta)
  )
}

#' @export
print.koine_family_spec <- function(x, ...) {
  cat("<koine family_spec>\n")
  cat("  name: ", x$name, "  link: ", x$link, "\n", sep = "")
  aux <- if (length(x$aux_params)) {
    paste(x$aux_params, collapse = ", ")
  } else {
    "(none)"
  }
  cat("  aux : ", aux, "\n", sep = "")
  invisible(x)
}
