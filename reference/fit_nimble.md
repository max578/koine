# Fit the NIMBLE-Gibbs lens

Runs NIMBLE's conjugate / Metropolis-within-Gibbs sampler on the pass's
GLM and returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md) of
equal-weight posterior draws with a complete provenance record. If the
family/link pair is unsupported, returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md)
rather than a wrong fit (I1).

## Usage

``` r
fit_nimble(
  pass,
  lens_spec = NULL,
  n_iter = 4000L,
  n_burn = 2000L,
  n_chains = 4L
)
```

## Arguments

- pass:

  A `koine_pass`.

- lens_spec:

  The lens's registry record; defaults to the registry's `nimble_gibbs`
  record.

- n_iter:

  Iterations per chain (including burn-in).

- n_burn:

  Burn-in iterations discarded per chain.

- n_chains:

  Number of chains.

## Value

A `koine_fb_fit`, or a `koine_abstention` for an unsupported link.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + 0.8 * d$x + rnorm(40)
pass <- koine_pass(y ~ x, d, gaussian())
fit_nimble(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
# }
```
