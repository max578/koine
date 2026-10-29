# Random-effects design for a pass (single or MET random intercepts)

Returns the per-observation group index/indices and group counts for the
grouping factor(s) of the random intercept(s). `NULL` for a population
pass. Used by the RE lenses (the marginal likelihood) and the integrated
predictive.

## Usage

``` r
.pass_re_design(pass)
```

## Arguments

- pass:

  A `koine_pass`.

## Value

A list as described above, or `NULL` when the pass has no random effect.

## Details

For the single random-intercept case the returned list carries the
legacy fields `grouping`, `group`, `levels`, and `n_groups` directly (so
the existing single-RE code path is byte-identical). It always also
carries `factors`, a per-grouping list of
`list(grouping, group, levels, n_groups)`, the `structure` (`"single"` /
`"crossed"` / `"nested"`), and the grouping names — the surface the MET
marginal core reads.
