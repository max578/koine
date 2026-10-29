# Fit the Polya-Gamma Gibbs lens (logistic regression)

Runs the exact Polya-Gamma Gibbs sampler for a Bernoulli-logit model and
returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md).
Resolves only for the binomial (logit) family; for any other family it
returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md)
rather than a wrong fit.

## Usage

``` r
fit_pg_gibbs(
  pass,
  lens_spec = NULL,
  n_iter = 3000L,
  n_burn = 1500L,
  n_chains = 4L
)
```

## Arguments

- pass:

  A `koine_pass` with a binomial (logit) family.

- lens_spec:

  The lens's registry record; defaults to the `pg_gibbs` record.

- n_iter:

  Iterations per chain.

- n_burn:

  Burn-in iterations discarded per chain.

- n_chains:

  Number of chains.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## References

Polson, Scott & Windle (2013). Bayesian inference for logistic models
using Polya-Gamma latent variables. *JASA* 108(504):1339-1349.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x = rnorm(80))
d$y <- rbinom(80, 1, plogis(0.5 * d$x))
pass <- koine_pass(y ~ x, d, binomial())
fit_pg_gibbs(pass, n_iter = 1000L, n_burn = 500L, n_chains = 2L)
# }
```
