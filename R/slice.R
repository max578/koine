## ----------------------------------------------------------------------------
## 1-D slice sampler --- the funnel-robust variance-component kernel
## (koine-authored; Neal 2003, "Slice Sampling", Ann. Statist.)
##
## A univariate slice sampler with stepping-out and shrinkage. koine uses it for
## ONE job: the variance-component coordinate (log tau) of the Poisson marginal
## mixed model, where the marginal posterior is a variance-component FUNNEL --- the
## conditioning of the (beta, log tau) geometry varies by position (the AGHQ-probed
## condition number runs 7 -> 124 across the log-tau range). A static-metric
## Hamiltonian kernel cannot serve the neck and the body at once (it collapses the
## step size, giving ~50% divergences and order-minutes fits). A slice sampler has
## no step size to mistune: it brackets the conditional automatically and so makes
## correct-scale moves at every point of the funnel. Blocked with the koine NUTS
## kernel on the coefficient sub-vector at FIXED tau (where the soft direction
## vanishes and the conditional is well-conditioned), it samples the marginal to
## calibration where the joint Hamiltonian kernel could not.
##
## It shares no code with any engine --- a self-contained, koine-authored slice
## kernel, used only as the variance-component block of the Hamiltonian-within-
## Gibbs Poisson scheme.
## ----------------------------------------------------------------------------

# One univariate slice-sampling transition (Neal 2003, Fig. 3 stepping-out +
# Fig. 5 shrinkage). Draws a new `x` from a density proportional to exp(logf(x))
# given the current `x0`, leaving that density invariant.
#
#   logf      the log of the (unnormalised) target density, x -> numeric
#   x0        current position (logf(x0) must be finite)
#   w         estimated typical slice width (the stepping-out unit)
#   m_step    cap on stepping-out expansions in each direction (bounds the work)
#   lower/upper  hard support bounds (default unbounded)
#
# Returns the new position. Each call uses the ambient RNG, so a set.seed upstream
# makes the whole chain reproducible (the determinism contract).
.slice_1d <- function(logf, x0, w = 1, m_step = 50L, lower = -Inf,
                      upper = Inf) {
  f0 <- logf(x0)
  if (!is.finite(f0)) {
    return(x0)                                  # degenerate; stay put (I2-safe)
  }
  # the vertical level defining the slice: log_y = f0 - Exp(1), so the slice is
  # { x : logf(x) > log_y }.
  log_y <- f0 - stats::rexp(1)

  # -- stepping out: grow an interval [l, r] around x0 until both ends are below
  #    the slice (or the step cap / support bound is hit).
  u <- stats::runif(1)
  l <- x0 - w * u
  r <- l + w
  j <- floor(m_step * stats::runif(1))
  k <- m_step - 1 - j
  while (j > 0 && l > lower && logf(l) > log_y) {
    l <- l - w
    j <- j - 1
  }
  while (k > 0 && r < upper && logf(r) > log_y) {
    r <- r + w
    k <- k - 1
  }
  l <- max(l, lower)
  r <- min(r, upper)

  # -- shrinkage: sample uniformly in [l, r], shrinking toward x0 on rejection
  #    until a point inside the slice is found. Guaranteed to terminate because
  #    the interval shrinks to x0, which is in the slice by construction.
  for (rep in seq_len(100L)) {
    x1 <- l + stats::runif(1) * (r - l)
    if (is.finite(logf(x1)) && logf(x1) > log_y) {
      return(x1)
    }
    if (x1 < x0) {
      l <- x1
    } else {
      r <- x1
    }
  }
  x0                                            # pathological; stay put
}
