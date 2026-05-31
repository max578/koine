# koine

> The synthesised fourth opinion for Bayesian triangulation.

## Description

**koine** synthesises one coherent, *independent* predictive opinion for a
Bayesian model from a curated battery of known and underutilised inference
methods --- or **abstains** when no diversity gain exists. It is designed to
serve as a fourth, synthesised witness alongside an independent panel of
inference engines.

The load-bearing idea: averaging *reduces* diversity and manufactures false
confidence. koine therefore synthesises on an axis **orthogonal** to the engine
axis --- predictive-space leave-one-out **stacking** (PSIS-gated) over
*model-structure and method* diversity, never averaging the engines. Its
disagreement carries a genuinely new signal ("sensitive to model structure /
inference family"), so it stays a real fourth witness rather than a weighted
echo of correlated ones.

koine is **programmatic infrastructure**: it exposes no user-facing modelling
verbs. A caller hands it an engine-neutral request and receives a synthesised
predictive opinion or a structured abstention.

## Installation

koine is in early development and is not yet distributed. From a local clone:

```r
# install.packages("pak")
pak::local_install(".")
```

## Usage

```r
library(koine)

# an engine-neutral request: formula, data, family, seed
d <- data.frame(y = rnorm(60), x = rnorm(60))
pass <- koine_pass(y ~ x, d, gaussian(), seed = 1L)

# synthesise the fourth opinion from the battery (or abstain)
res <- synthesise(pass)
res
```

`synthesise()` returns either a `koine_result` --- carrying the stacked
predictive opinion, the PSIS gate report, the diagnostic layer (weight
concentration, the over-dispersed linear-pool envelope, and the
prior-conditional evidence axis), and an independence report --- or a
first-class `koine_abstention`. Honest under-claiming beats a wrong fourth
opinion.

The battery spans three diversity axes: inference-family witnesses (DRAM,
NIMBLE-Gibbs, Polya-Gamma Gibbs, DREAM_ZS), a heavy-tail structure variant
(robust-t), and a shrinkage spectrum (ridge / lasso / horseshoe).

### Validation

koine-authored lenses must clear a correctness gate before they may corroborate:
simulation-based calibration and an exact grid-based reference-posterior check
(`inst/validation/`). An unvalidated lens abstains rather than vouch for itself.

## Dependencies

- **Imports** (always required): `digest`, `loo`, `posterior`, `RNetCDF`,
  `stats`, `utils`.
- **Suggests** (optional back-ends, gated at runtime): `BayesianTools`,
  `BayesLogit`, `bayesreg`, `bridgesampling`, `nimble`, `ncdf4` --- plus
  `knitr` / `rmarkdown` / `testthat` / `spelling` for the developer toolchain.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[Code of Conduct](CODE_OF_CONDUCT.md). In short: branch off `dev`, keep
`R CMD check --as-cran` clean, add a `NEWS.md` bullet for user-facing changes,
and reference the issue a pull request closes with `Fixes #<n>`.

## Licence

**No licence is declared yet.** The package is local-only and not distributed
while this is the case. A licence will be chosen only after the required
consultation among the interested parties has taken place; until then, please do
not redistribute koine in whole or in part. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff). Please also cite the methods koine drives ---
each lens's source paper is recorded in its documentation and provenance.
