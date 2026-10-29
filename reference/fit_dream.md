# Fit the DREAM_ZS lens

Runs BayesianTools' DREAMzs sampler on the pass's GLM and returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md).
Supports the Gaussian, Binomial (logit), and Poisson (log) families; for
any other family/link it returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md).

## Usage

``` r
fit_dream(pass, lens_spec = NULL, n_iter = 15000L, burn_frac = 0.5)
```

## Arguments

- pass:

  A `koine_pass`.

- lens_spec:

  The lens's registry record; defaults to the `dream_zs` record.

- n_iter:

  Total DREAM iterations (model evaluations across the ensemble).

- burn_frac:

  Fraction of each recorded chain discarded as burn-in.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## References

ter Braak & Vrugt (2008). Differential Evolution Markov Chain with
snooker updater and fewer chains. *Stat. Comput.* 18:435-446.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + 0.8 * d$x + rnorm(40)
pass <- koine_pass(y ~ x, d, gaussian())
fit_dream(pass, n_iter = 6000L)
# }
```
