# koine's opinion as a posterior draws array (the stacking mixture)

Materialises koine's single synthesised opinion as an equal-weight
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
over the parameters its contributing lenses share, by drawing from each
retained lens in proportion to its LOO-stacking weight. This is the
parameter-space representation a parent panel's parameter-marginal
triangulation consumes. The stacking weights, the stacked predictive
elpd, any per-lens parameters dropped from the shared set, and any scale
parameters whose lens set was restricted ride along as attributes
(`koine_weights`, `koine_elpd_stacked`, `koine_dropped_vars`,
`koine_scale_restricted`).

## Usage

``` r
koine_opinion_draws(result, n_draws = 4000L, seed = NULL)
```

## Arguments

- result:

  A `koine_result` produced with `synthesise(pass, keep_fits = TRUE)`. A
  `koine_abstention` is returned unchanged (koine contributed no opinion
  to represent).

- n_draws:

  Total draws in the returned mixture. Default 4000.

- seed:

  Optional integer seed for the mixture sampling; defaults to a seed
  derived from the pass's content hash, so the representation is
  reproducible.

## Value

A one-chain
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html),
or the input abstention.

## Details

The mixture is over koine's own battery (a legitimate single opinion),
never a pool across the parent's independent engines (charter SS4). Each
parameter is mixed only over the lenses for which it is comparable: the
linear-predictor coefficients (`b_<col>`) over every contributing lens,
but a scale parameter (`sigma`, `tau`) only over the family-likelihood
lenses — a likelihood-variant lens (the robust-t lens, which carries its
own Student-t likelihood) reports its scale on a different metric, so it
would bias the synthesised value if averaged in. The returned draws are
therefore a per-parameter marginal representation (each column is its
comparable-lens stacking mixture); the joint is approximate, which is
exactly what a parameter-marginal triangulation needs.

## See also

[`synthesise()`](https://max578.github.io/koine/reference/synthesise.md),
[`koine_can_contribute()`](https://max578.github.io/koine/reference/koine_can_contribute.md)

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = 1 + rnorm(80), x = rnorm(80))
pass <- koine_pass(y ~ x, d, gaussian())
res <- synthesise(pass, keep_fits = TRUE)
if (!is_abstention(res)) {
  dr <- koine_opinion_draws(res)
  posterior::summarise_draws(dr, "mean", "sd")
}
#> # A tibble: 3 × 3
#>   variable        mean     sd
#>   <chr>          <dbl>  <dbl>
#> 1 b_(Intercept)  1.08  0.103 
#> 2 b_x           -0.289 0.0952
#> 3 sigma          0.875 0.0667
# }
```
