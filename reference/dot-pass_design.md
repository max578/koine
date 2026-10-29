# Build the design matrix for a pass (fixed-effects terms)

Uses the pass's fixed-effects formula (random-effects terms stripped),
so the population-level design matrix is identical whether or not the
pass carries a random intercept.

## Usage

``` r
.pass_design(pass)
```

## Arguments

- pass:

  A `koine_pass`.

## Value

A list with `X` (model matrix), `y` (response vector), `N`, and `K`.
