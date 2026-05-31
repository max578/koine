# Contributing to koine

Thanks for your interest in koine. The package is in early development; the
notes below keep contributions consistent with the house conventions.

## Workflow

- Branch off `dev` (not `main`): `feature/<short-name>` or `bugfix/<short-name>`.
  `main` stays production-ready; merge upward only after the checks pass.
- Keep commits small and focused, and explain *why* in the message
  (e.g. "address #2 --- reconstruct pointwise log-likelihood via the family
  spec", not "minor change").
- A pull request that closes an issue must say so in its body: `Fixes #<n>`.

## Quality gates

Before opening a pull request:

```r
devtools::document()
devtools::test()                       # all green
rcmdcheck::rcmdcheck(args = "--as-cran")  # zero ERRORs / WARNINGs
```

- Add a [`NEWS.md`](NEWS.md) bullet for any user-facing change.
- New or changed code follows the two-space / 80-character house style; the
  repository carries an `air.toml` so [Air](https://posit-dev.github.io/air/)
  agrees with it. Do not restyle code unrelated to your change.
- Document every exported function with roxygen2 (markdown), including a runnable
  `@examples` block.

## Validating a new lens

A koine-authored lens must clear the correctness gate before it may contribute
to triangulation: simulation-based calibration and a reference-posterior check
(see `inst/validation/` and `koine_lens_validation_setup()`). Wire the lens's
generative model into the validation harness, run the study, record the verdict,
and attach the resulting `validation_ref` in the registry. An unvalidated lens
must abstain --- it is never silently trusted.

## Code of Conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
