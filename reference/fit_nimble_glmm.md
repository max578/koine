# Fit the NIMBLE GLMM lens

Runs NIMBLE's conjugate / Metropolis-within-Gibbs sampler on a mixed
model with a single random intercept and returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md) of
`(beta, sigma, tau)` draws — the random effects are sampled but not
retained, since the opinion is stacked on the integrated predictive. If
the family/link pair is unsupported, or the pass carries no supported
random intercept, returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md)
rather than a wrong fit (I1).

## Usage

``` r
fit_nimble_glmm(
  pass,
  lens_spec = NULL,
  n_iter = 4000L,
  n_burn = 2000L,
  n_chains = 4L
)
```

## Arguments

- pass:

  A `koine_pass` carrying a single random intercept `(1 | g)`, or the
  MET structure of two crossed/nested random intercepts
  `(1 | variety) + (1 | trial)`.

- lens_spec:

  The lens's registry record; defaults to the registry's `nimble_glmm`
  record.

- n_iter:

  Iterations per chain (including burn-in).

- n_burn:

  Burn-in iterations discarded per chain.

- n_chains:

  Number of chains.

## Value

A `koine_fb_fit`, or a `koine_abstention` for an unsupported pass.

## Examples

``` r
# \donttest{
set.seed(1)
g <- factor(rep(1:8, each = 5))
d <- data.frame(x = rnorm(40), g = g)
d$y <- 1 + 0.5 * d$x + rnorm(8)[g] + rnorm(40, 0, 0.7)
pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
fit_nimble_glmm(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
# }
```
