# koine

> The synthesised fourth opinion for Bayesian triangulation.

## Description

**koine** synthesises one coherent, *independent* predictive opinion for a
Bayesian model from a curated battery of known and underutilised inference
methods -- or **abstains** when no diversity gain exists. It is designed to
serve as a fourth, synthesised witness alongside an independent panel of
inference engines.

The load-bearing idea: averaging *reduces* diversity and manufactures false
confidence. koine therefore synthesises on an axis **orthogonal** to the engine
axis -- predictive-space leave-one-out **stacking** (PSIS-gated) over
*model-structure and method* diversity, never averaging the engines. Its
disagreement carries a genuinely new signal ("sensitive to model structure /
inference family"), so it stays a real fourth witness rather than a weighted
echo of correlated ones.

koine is **programmatic infrastructure**: it exposes no user-facing modelling
verbs. A caller hands it an engine-neutral request and receives a synthesised
predictive opinion or a structured abstention.

## Installation

koine is distributed through the private repository `max578/koine`:

```r
# install.packages("pak")
pak::pak("max578/koine")
```

or from a local clone:

```r
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

`synthesise()` returns either a `koine_result` -- carrying the stacked
predictive opinion, the PSIS gate report, the diagnostic layer (weight
concentration, the over-dispersed linear-pool envelope, and the
prior-conditional evidence axis), and an independence report -- or a
first-class `koine_abstention`. Stating less beats a wrong fourth
opinion.

The battery spans three diversity axes: inference-family witnesses (DRAM,
NIMBLE-Gibbs, Polya-Gamma Gibbs, DREAM_ZS), a heavy-tail structure variant
(robust-t), and a shrinkage spectrum (ridge / lasso / horseshoe).

koine also handles random *intercept* structure -- a single random intercept
`(1 | g)` and the breeder-standard multi-environment-trial (MET) structure of
two crossed or nested intercepts `(1 | variety) + (1 | trial)` -- on an axis
orthogonal to the sealed engines: NIMBLE-Gibbs *samples* the random effects while
the marginal lens *integrates* them out (exactly for the Gaussian LMM and the
Gaussian MET via a dense Cholesky of the marginal covariance, by adaptive
Gauss-Hermite quadrature for a single-intercept GLMM). The Gaussian, binomial,
and Poisson single-intercept GLMMs are enrolled; a non-Gaussian MET (a dense
crossed-GLMM marginal with no exact quadrature) abstains rather than approximate.
The MET lift lets koine corroborate a production MET model *directly* -- on the
variety-plus-trial structure that model fits -- rather than on a
single-random-effect reduction that discards the trial structure. A spectrum of
further lenses gives the random intercepts heavy
(Student-t) tails across a grid of degrees of freedom, so where the Gaussian-RE
lenses corroborate they earn stacking weight when group effects are outlying --
koine's RE opinion is sensitive to random-effect structure, not a weighted echo
of one assumption. All are stacked on the same per-observation integrated
predictive, which keeps the leave-one-out gate reliable; richer random-effects
structure (a random slope, more than two groupings) is refused
(`feature_unsupported`) rather than fitted as the wrong model.

### Validation

koine-authored lenses must clear a correctness gate before they may corroborate:
simulation-based calibration and an exact grid-based reference-posterior check
(`inst/validation/`). An unvalidated lens abstains rather than vouch for itself.

## Dependencies

- **Imports** (always required): `digest`, `loo`, `posterior`, `RNetCDF`,
  `stats`, `utils`.
- **Suggests** (optional back-ends, gated at runtime): `BayesianTools`,
  `BayesLogit`, `bayesreg`, `bridgesampling`, `nimble`, `ncdf4` -- plus
  `knitr` / `rmarkdown` / `testthat` / `spelling` for the developer toolchain.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[Code of Conduct](CODE_OF_CONDUCT.md). In short: branch off `main`, keep
`R CMD check --as-cran` clean, add a `NEWS.md` bullet for user-facing changes,
and reference the issue a pull request closes with `Fixes #<n>`.

## Licence

MIT. See [LICENSE](LICENSE) / [LICENSE.md](LICENSE.md).

## Citation

See [CITATION.cff](CITATION.cff). Please also cite the methods koine drives --
each lens's source paper is recorded in its documentation and provenance.
