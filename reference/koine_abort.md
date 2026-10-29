# Signal a koine contract-violation error

Raises a classed condition of subclass `koine_error_<code>` (and the
parent classes `koine_error` and `error`). The `code` must be one of the
closed error taxonomy in Stage-0 SS12; an unknown code is itself an
error, because a silent, untyped failure would violate charter hard-rule
2.

## Usage

``` r
koine_abort(code, message, ...)
```

## Arguments

- code:

  Single string; one of `.koine_error_codes`.

- message:

  Human-readable message.

- ...:

  Extra named fields attached to the condition for inspection.

## Value

Never returns; always stops.
