# Construct a structured abstention

An abstention records *why* koine declined to issue a fourth opinion, at
what scope, and which lenses (if any) were dropped on the way. It is
returned as a normal value and surfaces as "abstained" in the downstream
triangulation.

## Usage

``` r
abstention(
  reason,
  scope = c("synthesis", "lens"),
  detail,
  lenses_dropped = NULL
)
```

## Arguments

- reason:

  Single string from the closed taxonomy (Stage-0 SS8).

- scope:

  Either `"synthesis"` (the whole opinion) or `"lens"` (a single lens
  declined).

- detail:

  Human-readable explanation. Required — an abstention with no detail is
  itself misleading (I1).

- lenses_dropped:

  Optional list of `list(id =, reason =)` records for lenses dropped
  before the abstention.

## Value

An object of class `koine_abstention`.

## Examples

``` r
abstention("no_diverse_battery", detail = "only 1 opinion lens resolved")
#> <koine abstention>
#>   scope : synthesis
#>   reason: no_diverse_battery
#>   detail: only 1 opinion lens resolved
```
