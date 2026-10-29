# Fit the marginal mixed-model lens

Runs the koine-authored No-U-Turn Sampler on the random-effects
*marginal* model — the random intercepts integrated out, exactly for the
Gaussian LMM and by adaptive Gauss-Hermite quadrature for a GLMM — and
returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md) of
`(beta, sigma, tau)` draws. For the Poisson marginal — whose
variance-component funnel a single static-metric kernel cannot traverse
— a funnel-safe blocked sampler is used: `log tau` is slice-sampled, and
the coefficients are slice-sampled along their `tau`-resolved principal
axes (a Gibbs block that needs no step size and cannot stall); the
Gaussian and binomial marginals keep the joint NUTS kernel.

## Usage

``` r
fit_glmm_marginal(
  pass,
  lens_spec = NULL,
  n_iter = 1500L,
  n_burn = 750L,
  n_chains = 4L,
  n_nodes = 20L,
  max_treedepth = 9L
)
```

## Arguments

- pass:

  A `koine_pass` carrying a single random intercept `(1 | g)`, or the
  MET structure of two crossed/nested random intercepts
  `(1 | variety) + (1 | trial)`.

- lens_spec:

  The lens's registry record; defaults to the registry's `glmm_marginal`
  record.

- n_iter:

  Iterations per chain (including warm-up). The Hamiltonian kernel needs
  far fewer than a random walk for the same effective sample size, so
  the default is modest.

- n_burn:

  Warm-up iterations per chain (step-size + mass-matrix adaptation),
  discarded from the returned draws.

- n_chains:

  Number of chains.

- n_nodes:

  Gauss-Hermite nodes for the GLMM marginal (ignored for the exact
  Gaussian LMM path).

- max_treedepth:

  Maximum NUTS tree depth (the joint-NUTS Gaussian / binomial path
  only), capping the leapfrog budget per iteration at `2^max_treedepth`;
  saturation of the cap is surfaced in the provenance trace, never
  silently absorbed. The Poisson blocked slice kernel has no tree depth.

## Value

A `koine_fb_fit`, or a `koine_abstention` for an unsupported pass.

## Details

For the breeder-standard multi-environment-trial (MET) structure — two
crossed or nested random intercepts `(1 | variety) + (1 | trial)` — the
Gaussian marginal stays exact (a dense Cholesky of
`V = sigma^2 I + sum_r tau_r^2 Z_r Z_r'`) and is sampled by a
coordinate-wise slice kernel over
`(beta, sigma, tau_variety, tau_trial)`; the draws then carry one
`tau_<grouping>` per random factor. A non-Gaussian MET (a dense
crossed-GLMM marginal, which has no exact quadrature) returns a
lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md).
The pointwise log-likelihood is left `NULL` for
[`assemble_loglik()`](https://max578.github.io/koine/reference/assemble_loglik.md)
to compute as the per-observation integrated predictive (F3). If the
pass carries no supported random intercept, returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md)
rather than a wrong fit (I1).

## Examples

``` r
# \donttest{
set.seed(1)
g <- factor(rep(1:8, each = 5))
d <- data.frame(x = rnorm(40), g = g)
d$y <- 1 + 0.5 * d$x + rnorm(8)[g] + rnorm(40, 0, 0.7)
pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
fit_glmm_marginal(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
# }
```
