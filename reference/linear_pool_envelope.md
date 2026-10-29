# Linear-pool envelope diagnostic (Stage-0 SS9.2)

An equal-weight mixture of the contributing lenses' conditional-mean
predictives, deliberately **over-dispersed** to serve as a reference
bracket rather than a confident prediction. For each observation the
envelope half-width combines the average within-lens posterior
uncertainty with the between-lens disagreement, the latter inflated by
`inflate` (default 2). The reported envelope is the full equal-weight
pool; flagging is leave-one-out — a lens is flagged at an observation
when its posterior-mean response falls outside the over-dispersed
bracket formed by the *other* lenses, so a disagreement carrier is
caught against the consensus rather than inflating its own bracket. The
envelope is a diagnostic surface only and is never shipped as the
opinion (I2).

## Usage

``` r
linear_pool_envelope(fb_fits, pass, inflate = 2, level = 0.95)
```

## Arguments

- fb_fits:

  A list of contributed `koine_fb_fit` objects (opinion lenses).

- pass:

  The originating `koine_pass`.

- inflate:

  Over-dispersion multiplier on the between-lens disagreement variance
  (default 2, per SS9.2's "inflated approximately 2x disagreement").

- level:

  Central interval width for the bracket (default 0.95).

## Value

A list with `obs_id`, `pool_mean`, `lower`, `upper` (length-N vectors),
`flags` (a data frame of `lens_id` x out-of-envelope counts),
`n_components`, and `inflate`; or `NULL` when fewer than two lenses
expose a comparable mean-response (no envelope is meaningful).

## Details

Operates on the conditional-mean scale `mu_i = linkinv(eta_i)`, which is
comparable across every lens regardless of likelihood-tail differences,
so a robust-t lens and a Gaussian lens are bracketed on the same
footing.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = 1 + rnorm(50), x = rnorm(50))
pass <- koine_pass(y ~ x, d, gaussian())
fits <- list(
  fit_dram(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L),
  fit_robust_t(pass, n_iter = 800L, n_burn = 400L, n_chains = 2L))
env <- linear_pool_envelope(fits, pass)
env$flags
#>    lens_id n_outside n_obs
#> 1     dram         0    50
#> 2 robust_t         0    50
# }
```
