# Per-observation integrated predictive for a crossed/nested Gaussian MET model

For a *new* observation the random effect of every grouping is drawn
afresh and marginalised, so the predictive is
`N(y_i | eta_fixed_i, sigma^2 + sum_r tau_r^2)` — the same held-out
target every MET lens is stacked on, keeping the leave-one-out PSIS gate
in its reliable regime (the analogue of the single-RE integrated
predictive).

## Usage

``` r
.re_integrated_loglik_met(y, eta_fixed, taus, sigma)
```

## Arguments

- y:

  Response vector (length `N`).

- eta_fixed:

  Fixed-part linear predictor (length `N`).

- taus:

  Vector of random-intercept SDs (one per grouping).

- sigma:

  Residual SD.

## Value

A length-`N` vector of integrated pointwise log-likelihoods.
