# Is this object an abstention?

Is this object an abstention?

## Usage

``` r
is_abstention(x)
```

## Arguments

- x:

  Any object.

## Value

`TRUE` if `x` is a `koine_abstention`.

## Examples

``` r
out <- abstention("no_diverse_battery", scope = "synthesis",
                  detail = "only one opinion lens resolved.")
is_abstention(out)
#> [1] TRUE
is_abstention(42)
#> [1] FALSE
```
