# Resolve the battery for a pass (family-conditional enrolment, Stage-0 SS3.3)

Enrols a lens iff: its tier is allowed; its family predicate passes; its
feature support covers the formula; the log-likelihood is producible
(native or via the family spec); the adapter is implemented in this
build; and, for a koine-authored lens, its validation evidence is on
record (the open-arms correctness gate, SS11 — an unvalidated
first-mover abstains). A lens that fails any test is excluded *with a
recorded reason* — a logged absence, not an error (I1).

## Usage

``` r
resolve_battery(pass, registry = koine_stage1_registry())
```

## Arguments

- pass:

  A `koine_pass`.

- registry:

  A registry list; defaults to
  [`koine_stage1_registry()`](https://max578.github.io/koine/reference/koine_stage1_registry.md).

## Value

A list with `enrolled` (lens specs) and `excluded`
(`list(id =, reason =)` records).

## Examples

``` r
set.seed(1)
d <- data.frame(y = rnorm(30), x = rnorm(30))
pass <- koine_pass(y ~ x, d, gaussian())
resolved <- resolve_battery(pass)
vapply(resolved$enrolled, function(L) L$id, character(1))
#> [1] "dram"         "nimble_gibbs" "dream_zs"     "robust_t"     "horseshoe"   
#> [6] "lasso"        "ridge"       
```
