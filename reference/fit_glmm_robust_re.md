# Fit the heavy-tailed random-effect lens (Student-t random intercepts)

Fits a Gaussian linear mixed model whose random intercepts are Student-t
(`b_g ~ t_nu(0, tau)`) by delayed-rejection adaptive Metropolis on the
t-RE marginal, and returns an
[`fb_fit()`](https://max578.github.io/koine/reference/fb_fit.md)
carrying its own t-RE integrated predictive. Resolves only for a
Gaussian (identity) response with a supported single random intercept;
otherwise returns a lens-scope
[`abstention()`](https://max578.github.io/koine/reference/abstention.md).

## Usage

``` r
fit_glmm_robust_re(
  pass,
  lens_spec = NULL,
  nu = 4,
  n_iter = 4000L,
  n_burn = 2000L,
  n_chains = 4L,
  n_nodes = 32L
)
```

## Arguments

- pass:

  A `koine_pass` with a Gaussian family and a single `(1 | g)`.

- lens_spec:

  The lens's registry record; defaults to the `glmm_robust_re` record.

- nu:

  Degrees of freedom of the Student-t random effect (fixed; lower =
  heavier-tailed group effects). Default 4.

- n_iter, n_burn, n_chains:

  Sampler controls.

- n_nodes:

  Gauss-Laguerre nodes for the scale-mixture integral.

## Value

A `koine_fb_fit`, or a `koine_abstention`.

## Examples

``` r
# \donttest{
set.seed(1)
g <- factor(rep(1:12, each = 8))
bt <- stats::rt(12, df = 4) * 1.2
d <- data.frame(x = rnorm(96), g = g)
d$y <- 1 + 0.5 * d$x + bt[g] + rnorm(96, 0, 0.7)
pass <- koine_pass(y ~ x + (1 | g), d, gaussian())
fit_glmm_robust_re(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L)
# }
```
