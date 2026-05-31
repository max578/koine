# Write a grouped NetCDF-4 fixture with RNetCDF (R-side, no Python). The
# cross-LANGUAGE direction (Python write -> R read) is covered durably by
# spikes/netcdf_roundtrip/; here we unit-test the reader's group traversal,
# name-based dim resolution (A2), and the hard assertions A1-A4.
write_fixture_nc <- function(path, n_chain = 2L, n_draw = 10L, n_obs = 4L,
                             n_beta = 2L, obs_id = NULL) {
  if (is.null(obs_id)) obs_id <- sprintf("obs_%04d", seq_len(n_obs))
  nc <- RNetCDF::create.nc(path, format = "netcdf4")
  on.exit(RNetCDF::close.nc(nc), add = TRUE)

  prov <- list(
    method = "fixture_lens", method_version = "0.0.1", impl_language = "Python",
    seed = 11L, derived_seed = 22L, ir_hash = "fixturehash",
    timestamp_iso8601 = "2026-05-31T00:00:00+0000",
    maturity_tier = "reference-validated", role = "opinion",
    exactness = "exact", independence_family = "nested",
    independence_kernel = "dynesty")
  put_prov <- function(grp) {
    for (nm in names(prov)) {
      v <- prov[[nm]]
      if (is.numeric(v)) {
        RNetCDF::att.put.nc(grp, "NC_GLOBAL", nm, "NC_INT", as.integer(v))
      } else {
        RNetCDF::att.put.nc(grp, "NC_GLOBAL", nm, "NC_CHAR", v)
      }
    }
  }
  put_prov(nc)                                   # root mirror (for ncdf4 path)

  # -- posterior group --
  gp <- RNetCDF::grp.def.nc(nc, "posterior")
  RNetCDF::dim.def.nc(gp, "chain", n_chain)
  RNetCDF::dim.def.nc(gp, "draw", n_draw)
  RNetCDF::dim.def.nc(gp, "beta_dim", n_beta)
  RNetCDF::var.def.nc(gp, "beta", "NC_DOUBLE", c("chain", "draw", "beta_dim"))
  RNetCDF::var.def.nc(gp, "sigma", "NC_DOUBLE", c("chain", "draw"))
  beta <- array(rnorm(n_chain * n_draw * n_beta),
                dim = c(n_chain, n_draw, n_beta))
  RNetCDF::var.put.nc(gp, "beta", beta)
  RNetCDF::var.put.nc(gp, "sigma",
                      array(abs(rnorm(n_chain * n_draw)) + 0.5,
                            dim = c(n_chain, n_draw)))
  put_prov(gp)

  # -- log_likelihood group: SCRAMBLED dim order (draw, obs_id, chain) --
  gl <- RNetCDF::grp.def.nc(nc, "log_likelihood")
  RNetCDF::dim.def.nc(gl, "chain", n_chain)
  RNetCDF::dim.def.nc(gl, "draw", n_draw)
  RNetCDF::dim.def.nc(gl, "obs_id", n_obs)
  RNetCDF::var.def.nc(gl, "obs", "NC_DOUBLE", c("draw", "obs_id", "chain"))
  RNetCDF::var.def.nc(gl, "obs_id", "NC_STRING", "obs_id")
  # canonical[c,d,o] = c*1000 + d*100 + o ; store as (draw, obs_id, chain)
  canon <- array(0, dim = c(n_chain, n_draw, n_obs))
  for (c in seq_len(n_chain)) for (d in seq_len(n_draw)) for (o in seq_len(n_obs))
    canon[c, d, o] <- c * 1000 + d * 100 + o
  RNetCDF::var.put.nc(gl, "obs", aperm(canon, c(2, 3, 1)))
  RNetCDF::var.put.nc(gl, "obs_id", obs_id)

  # -- observed_data group --
  go <- RNetCDF::grp.def.nc(nc, "observed_data")
  RNetCDF::dim.def.nc(go, "obs_id", n_obs)
  RNetCDF::var.def.nc(go, "y", "NC_DOUBLE", "obs_id")
  RNetCDF::var.def.nc(go, "obs_id", "NC_STRING", "obs_id")
  RNetCDF::var.put.nc(go, "y", rnorm(n_obs))
  RNetCDF::var.put.nc(go, "obs_id", obs_id)
  path
}

test_that("read_idata_nc (RNetCDF) traverses groups and resolves A1-A4", {
  skip_if_not_installed("RNetCDF")
  path <- tempfile(fileext = ".nc")
  obs_id <- sprintf("obs_%04d", 1:4)
  write_fixture_nc(path, obs_id = obs_id)

  fit <- read_idata_nc(path, reader = "rnetcdf")
  expect_s3_class(fit, "koine_fb_fit")
  expect_identical(fit$status, "contributed")

  # A1 + A2: log_likelihood is (chain, draw, obs_id) after name-based reorder.
  ll <- fit$log_likelihood
  expect_equal(dim(ll), c(2L, 10L, 4L))
  expect_equal(ll[1L, 1L, 1L], 1 * 1000 + 1 * 100 + 1)
  expect_equal(ll[2L, 10L, 4L], 2 * 1000 + 10 * 100 + 4)

  # posterior draws flattened to indexed variables.
  expect_setequal(posterior::variables(fit$draws),
                  c("beta[1]", "beta[2]", "sigma"))

  # provenance read from group attrs.
  expect_identical(fit$provenance$method, "fixture_lens")
  expect_identical(fit$provenance$independence_key$family, "nested")
})

test_that("read_idata_nc enforces the obs-axis assertion (A3)", {
  skip_if_not_installed("RNetCDF")
  path <- tempfile(fileext = ".nc")
  write_fixture_nc(path, obs_id = sprintf("obs_%04d", 1:4))
  d <- data.frame(y = rnorm(4), x = rnorm(4))
  pass_ok <- koine_pass(y ~ x, d, gaussian(), obs_id = sprintf("obs_%04d", 1:4))
  expect_s3_class(read_idata_nc(path, pass = pass_ok, reader = "rnetcdf"),
                  "koine_fb_fit")

  pass_bad <- koine_pass(y ~ x, d, gaussian(), obs_id = c("z1", "z2", "z3", "z4"))
  expect_error(read_idata_nc(path, pass = pass_bad, reader = "rnetcdf"),
               class = "koine_error_obs_axis_mismatch")
})

test_that("ncdf4 fallback reads the same file via root-mirrored provenance", {
  skip_if_not_installed("RNetCDF")
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  write_fixture_nc(path, obs_id = sprintf("obs_%04d", 1:4))

  fit <- read_idata_nc(path, reader = "ncdf4")
  expect_s3_class(fit, "koine_fb_fit")
  expect_equal(dim(fit$log_likelihood), c(2L, 10L, 4L))
  expect_equal(fit$log_likelihood[2L, 10L, 4L], 2 * 1000 + 10 * 100 + 4)
  expect_identical(fit$provenance$method, "fixture_lens")
})
