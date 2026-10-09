# Integrating koine as a fourth inference backend

This note specifies the contract a parent triangulation package (flexyBayes)
drives to use **koine** as a fourth, synthesised witness alongside an
independent panel of inference engines. koine is human-directly-unusable
infrastructure: the parent hands it a model and a required output format, and
koine returns one synthesised opinion or a structured abstention. The parent
owns the IR translation and the backend wrapper; koine stays
parent-agnostic (it has no dependency on, and no knowledge of, the parent's
classes).

## What koine contributes

A single, synthesised, independent predictive opinion formed by
predictive-space LOO-stacking over a curated battery of inference lenses
(orthogonal to the engine axis), or a first-class abstention when no diversity
gain exists. koine is one opinion versus the parent's independent natives; it is
never triangulated against itself, and its stacking-weights-over-lenses
diagnostic is surfaced, never shipped as a pooled posterior over the engines.

## The supported-model corner (capability)

koine contributes only where its battery can field a diverse, reliable opinion.
In this build (0.3.0) that is: the population GLM; the single random intercept
`(1 | g)` on one grouping factor; and the **MET structure** --- two crossed or
nested random intercepts (`(1 | variety) + (1 | trial)`) --- for the
`gaussian`, `binomial`, and `poisson` families (Poisson enrolled at 0.2.0,
validated by SBC + NIMBLE concordance). For the Gaussian MET, the marginal is
exact (dense Cholesky); a non-Gaussian MET fits via the marginalised
NUTS/blocked-slice kernel and can itself abstain (`convergence_failed`) where
the marginalised variance component is too weakly identified to mix, per
lens (`R/registry.R`, the `glmm_marginal` record).

Everything richer --- random slopes, three or more crossed groupings,
correlated effects, factor-analytic or genomic (`vm` / `ped`) structures,
splines, `AR1`, GxE, multivariate responses --- is **out of corner**, and
koine abstains rather than mis-modelling it. Abstention is the expected default,
not a failure.

Use the pre-flight predicate to gate availability cheaply (it resolves the
battery but fits nothing):

```r
ok <- koine::koine_can_contribute(formula, data, family)
# TRUE, or a reason string: "feature_unsupported" / "no_diverse_battery" /
# "pass_invalid". A TRUE is the minimum condition (the PSIS gate and weight-
# degeneracy check are only knowable after fitting, so koine may still abstain
# at synthesis time).
```

## The handshake

1. **Translate the parent IR to a koine request.** Use the IR-targetable builder
   `koine_pass_from_design()` --- it takes a named design matrix, a response, a
   family, and an optional grouping vector directly (no formula, re-entrant,
   deparse-free), so the parent hands koine its own design arrays:

   | parent IR (`fb_terms` / design)            | koine request (`koine_pass_from_design`) |
   |--------------------------------------------|------------------------------------------|
   | the fixed-effects design matrix (named cols, incl. `(Intercept)`) | `x` |
   | the response                               | `y`                                      |
   | `family`, `link`                           | `family`                                 |
   | a single random intercept's grouping factor | `group`                                 |
   | two crossed/nested random intercepts (MET) | `groups = list(variety = v, trial = t)` |
   | any out-of-corner term                     | koine abstains (or pre-gate via `koine_can_contribute()`) |

   ```r
   pass <- koine::koine_pass_from_design(x, y, family, group = g,
                                         seed = <integer>)
   # or, for the MET structure (supply `group` OR `groups`, never both):
   pass <- koine::koine_pass_from_design(x, y, family,
                                         groups = list(variety = v, trial = t),
                                         seed = <integer>)
   ```

   (A formula entry, `koine::koine_pass(formula, data, family, ...)`, is also
   available for human / exploratory use.)

2. **Synthesise the opinion, retaining the fits for the draws representation:**

   ```r
   res <- koine::synthesise(pass, keep_fits = TRUE)
   ```

   `res` is a `koine_result` or a `koine_abstention` --- never a silent partial.

3. **Branch on the outcome.**

   ```r
   if (koine::is_abstention(res)) {
     # res$reason (closed taxonomy), res$detail, res$lenses_dropped
     # -> the parent records koine's column as "abstained"
   } else {
     # res$opinion$weights        stacking weights over the battery (a simplex)
     # res$opinion$elpd_stacked   the stacked LOO predictive density
     # res$diagnostics            weight_concentration / linear_pool_envelope /
     #                            evidence_axis (surfaced, never pooled)
     # res$independence           the per-lens method + diversity-axis report
   }
   ```

4. **Materialise the parameter-space representation** for a parameter-marginal
   triangulation (the parent's `fb_as_draws_simple`):

   ```r
   dr <- koine::koine_opinion_draws(res)   # a posterior::draws_array, one chain
   ```

   `dr` is koine's opinion as draws over the parameters its contributing lenses
   share, sampled from each lens in proportion to its stacking weight. It carries
   the predictive-space signal as attributes so the parent never loses it:
   `attr(dr, "koine_weights")`, `attr(dr, "koine_elpd_stacked")`,
   `attr(dr, "koine_dropped_vars")`. A `koine_abstention` passes through
   unchanged.

   **Caveat (load-bearing).** The mixture is over koine's own battery --- a
   legitimate representation of koine's single opinion --- and is **not** a pool
   across the parent's independent engines (that would manufacture false
   confidence). Each parameter is mixed only over the lenses for which it is
   comparable: the fixed-effect coefficients over every contributing lens, but a
   scale parameter (`sigma`, `tau`) only over the family-likelihood lenses --- a
   likelihood-variant lens (the robust-t lens) reports its scale on a different
   metric, so it is excluded from the scale mixture (and flagged in the
   `koine_scale_restricted` attribute) while still contributing to the
   coefficients. The draws are thus a per-parameter marginal representation; the
   joint is approximate, which is what a parameter-marginal triangulation needs.

## Parameter names (for the canonical-name mapper)

koine's draws use stable variable names the parent maps to its canonical
(brms-style) scheme:

| koine variable     | meaning                          | canonical (brms-style) |
|--------------------|----------------------------------|------------------------|
| `b_(Intercept)`    | fixed intercept                  | `b_Intercept`          |
| `b_<col>`          | fixed coefficient for `<col>`    | `b_<col>`              |
| `sigma`            | Gaussian residual SD             | `sigma`                |
| `tau`              | single random-intercept SD       | `sd_<g>__Intercept`    |
| `tau_<grouping>`   | MET random-intercept SD (one per crossed/nested factor) | `sd_<grouping>__Intercept` |

## Independence registration (parent ADR 0029)

koine is not a single inference paradigm --- it is a *stacked battery* spanning
several (random-walk, Gibbs, Hamiltonian, population, tempering, ...). Two
consequences for the parent's independence registry:

* Register koine's pairs with each engine on the **algorithmic** and
  **implementation** axes: koine shares no code with any engine, and its
  synthesised opinion is not a single mainstream algorithm.
* koine's *realised* independence is weight-dependent: when its stacking weight
  rests on a lens whose family coincides with an engine (e.g. the
  `hamiltonian` marginal lens versus an HMC engine), the algorithmic-axis
  independence weakens for that opinion, though implementation independence
  holds (koine's NUTS shares no code with Stan). The parent should read
  `res$independence` to surface koine's realised, weight-weighted independence
  profile rather than a static label.

## Status and governance gates

* **Version** `0.3.0` (released; private `max578/koine`). The parent gates
  availability on koine's presence (a `Suggests` dependency).
* **Licence** `MIT + file LICENSE` (declared; see `LICENSE`/`LICENSE.md`).
* **Cross-engine concordance** (the third validation criterion) is deferred to
  graduation: koine's lenses must be concordance-checked against the engine
  panel on shared in-corner models once the parent is present.
* The founding orthogonal-diversity principle is recorded for a blinded
  consensus to ratify at graduation.
