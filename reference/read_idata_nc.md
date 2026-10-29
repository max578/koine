# Read an ArviZ DataTree NetCDF-4 file into an `fb_fit`

Read an ArviZ DataTree NetCDF-4 file into an `fb_fit`

## Usage

``` r
read_idata_nc(path, pass = NULL, reader = c("auto", "rnetcdf", "ncdf4"))
```

## Arguments

- path:

  Path to the `.nc` file.

- pass:

  Optional `koine_pass`; when supplied, the log-likelihood observation
  axis is asserted equal to `pass$obs_dim_names$obs_id` in the same
  order (A3).

- reader:

  One of `"auto"` (RNetCDF if available, else ncdf4), `"rnetcdf"`, or
  `"ncdf4"`.

## Value

A `koine_fb_fit`.

## Examples

``` r
# read_idata_nc() resolves an ArviZ DataTree NetCDF written by a Tier-3
# (Python) lens into the uniform fb_fit contract; supply the .nc path and
# the originating pass:
if (FALSE) { # \dontrun{
fit <- read_idata_nc("inference.nc", pass)
} # }
```
