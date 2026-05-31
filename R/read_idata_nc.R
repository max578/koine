## ----------------------------------------------------------------------------
## read_idata_nc --- the dual-reader NetCDF interchange (Stage-0 SS5)
##
## Reads an ArviZ-style DataTree NetCDF-4 file (one group per fb_fit field) into
## an `fb_fit`, with no live Python at read time. Two backends, chosen
## empirically (spike `spikes/netcdf_roundtrip/`):
##   * canonical: RNetCDF (reads groups + group-level attrs via the group API).
##   * fallback : ncdf4   (reads grouped data via `group/var` path-names, root
##                global attrs only --- relies on provenance being mirrored to
##                root).
## Both resolve dims by NAME (A2) and enforce the four hard assertions A1-A4.
## ----------------------------------------------------------------------------

.bare <- function(x) sub("^.*/", "", x)

# Provenance attribute names the reader reconstructs from group/root attrs.
.prov_attr_names <- c(
  "method", "method_version", "impl_language", "seed", "derived_seed",
  "ir_hash", "timestamp_iso8601", "maturity_tier", "role", "exactness",
  "independence_family", "independence_kernel"
)

#' Read an ArviZ DataTree NetCDF-4 file into an `fb_fit`
#'
#' @param path Path to the `.nc` file.
#' @param pass Optional `koine_pass`; when supplied, the log-likelihood
#'   observation axis is asserted equal to `pass$obs_dim_names$obs_id` in the
#'   same order (A3).
#' @param reader One of `"auto"` (RNetCDF if available, else ncdf4),
#'   `"rnetcdf"`, or `"ncdf4"`.
#' @return A `koine_fb_fit`.
#' @examples
#' # read_idata_nc() resolves an ArviZ DataTree NetCDF written by a Tier-3
#' # (Python) lens into the uniform fb_fit contract; supply the .nc path and
#' # the originating pass:
#' \dontrun{
#' fit <- read_idata_nc("inference.nc", pass)
#' }
#' @export
read_idata_nc <- function(path, pass = NULL,
                          reader = c("auto", "rnetcdf", "ncdf4")) {
  reader <- match.arg(reader)
  stopifnot(file.exists(path))
  if (reader == "auto") {
    reader <- if (requireNamespace("RNetCDF", quietly = TRUE)) "rnetcdf"
              else "ncdf4"
  }
  groups <- switch(reader,
                   rnetcdf = .nc_read_rnetcdf(path),
                   ncdf4 = .nc_read_ncdf4(path))
  .assemble_fb_fit_from_groups(groups, pass)
}

# -- RNetCDF backend ---------------------------------------------------------
.nc_read_rnetcdf <- function(path) {
  nc <- RNetCDF::open.nc(path)
  on.exit(RNetCDF::close.nc(nc), add = TRUE)
  root <- RNetCDF::grp.inq.nc(nc)
  child_ids <- root$grps
  gnames <- vapply(child_ids, function(g) RNetCDF::grp.inq.nc(g)$name,
                   character(1))
  groups <- list()
  for (i in seq_along(child_ids)) {
    g <- child_ids[[i]]
    ginfo <- RNetCDF::grp.inq.nc(g)
    vars <- list()
    for (vid in ginfo$varids) {
      vinq <- RNetCDF::var.inq.nc(g, vid)
      dnames <- vapply(vinq$dimids,
                       function(d) RNetCDF::dim.inq.nc(g, d)$name, character(1))
      vars[[vinq$name]] <- list(
        data = RNetCDF::var.get.nc(g, vinq$name),
        dims = .bare(dnames),
        type = vinq$type
      )
    }
    groups[[gnames[i]]] <- list(
      vars = vars,
      attrs = .read_attrs_rnetcdf(g)
    )
  }
  groups
}

.read_attrs_rnetcdf <- function(grpid) {
  out <- list()
  for (nm in .prov_attr_names) {
    out[[nm]] <- tryCatch(
      RNetCDF::att.get.nc(grpid, "NC_GLOBAL", nm),
      error = function(e) NULL)
  }
  out
}

# -- ncdf4 backend (fallback) ------------------------------------------------
.nc_read_ncdf4 <- function(path) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("read_idata_nc(reader = 'ncdf4'): the ncdf4 package is not installed.",
         call. = FALSE)
  }
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  paths <- names(nc$var)
  gnames <- unique(sub("/[^/]*$", "", paths))
  root_attrs <- ncdf4::ncatt_get(nc, 0)
  groups <- list()
  for (gn in gnames) {
    members <- paths[sub("/[^/]*$", "", paths) == gn]
    vars <- list()
    for (p in members) {
      v <- nc$var[[p]]
      vars[[.bare(p)]] <- list(
        data = ncdf4::ncvar_get(nc, p),
        dims = .bare(vapply(v$dim, function(d) d$name, character(1))),
        type = v$prec
      )
    }
    # ncdf4 reads only root global attrs (the provenance-mirroring rule, SS5).
    attrs <- root_attrs[intersect(names(root_attrs), .prov_attr_names)]
    groups[[gn]] <- list(vars = vars, attrs = attrs)
  }
  groups
}

# -- shared assembler --------------------------------------------------------
.assemble_fb_fit_from_groups <- function(groups, pass) {
  if (is.null(groups[["posterior"]])) {
    koine_abort("netcdf_shape", "no 'posterior' group in the NetCDF file.")
  }
  draws <- .group_to_draws(groups[["posterior"]], require_chain_draw = TRUE)

  log_likelihood <- NULL
  if (!is.null(groups[["log_likelihood"]])) {
    ll_grp <- groups[["log_likelihood"]]
    vnm <- names(ll_grp$vars)[1L]
    ll <- .reorder_by_name(ll_grp$vars[[vnm]],
                           c("chain", "draw", "obs_id"))
    if (!identical(ll$type, "NC_DOUBLE") && !identical(ll$type, "double")) {
      koine_abort("netcdf_dtype",
                  sprintf("log_likelihood dtype is %s, expected float64.",
                          ll$type))
    }
    if (!is.null(pass)) {
      obs_read <- as.character(ll_grp$vars[["obs_id"]]$data)
      if (!identical(obs_read, pass$obs_dim_names$obs_id)) {
        koine_abort("obs_axis_mismatch",
                    "log_likelihood obs axis != pass$obs_dim_names$obs_id.")
      }
    }
    log_likelihood <- ll$data
  }

  sample_stats <- NULL
  if (!is.null(groups[["sample_stats"]])) {
    sample_stats <- .group_to_draws(groups[["sample_stats"]],
                                    require_chain_draw = TRUE)
  }

  prov <- .attrs_to_provenance(groups[["posterior"]]$attrs)
  fb_fit(draws = draws, provenance = prov, exactness = prov$exactness,
         sample_stats = sample_stats, log_likelihood = log_likelihood,
         status = "contributed")
}

# Reorder a variable's array so its named dims match `want` (A2), strictly by
# name. Errors with netcdf_dimorder if a wanted dim is absent.
.reorder_by_name <- function(var, want) {
  if (!setequal(var$dims, want)) {
    koine_abort("netcdf_dimorder",
                sprintf("dims {%s} do not match expected {%s}.",
                        paste(var$dims, collapse = ","),
                        paste(want, collapse = ",")))
  }
  perm <- match(want, var$dims)
  list(data = aperm(var$data, perm), type = var$type)
}

# Convert a posterior/sample_stats group into a posterior draws_array,
# flattening any extra parameter dims into indexed variable names.
.group_to_draws <- function(group, require_chain_draw = TRUE) {
  var_arrays <- list()
  n_chain <- NA_integer_
  n_draw <- NA_integer_
  for (vn in names(group$vars)) {
    v <- group$vars[[vn]]
    if (!all(c("chain", "draw") %in% v$dims)) {
      if (require_chain_draw) next
    }
    extra <- setdiff(v$dims, c("chain", "draw"))
    a <- .reorder_by_name(v, c("chain", "draw", extra))$data
    dd <- dim(a)
    if (is.na(n_chain)) {
      n_chain <- dd[1L]
      n_draw <- dd[2L]
    }
    if (length(extra) == 0L) {
      var_arrays[[vn]] <- array(a, dim = c(dd[1L], dd[2L]))
    } else {
      idx <- expand.grid(lapply(dd[-(1:2)], seq_len))
      for (r in seq_len(nrow(idx))) {
        lab <- paste0(vn, "[", paste(unlist(idx[r, ]), collapse = ","), "]")
        sub <- a
        # index the trailing dims
        args <- c(list(sub), list(TRUE, TRUE), as.list(unlist(idx[r, ])))
        var_arrays[[lab]] <- do.call(`[`, args)
      }
    }
  }
  # A1: nchain * ndraw == total draws per variable
  vn1 <- var_arrays[[1L]]
  if (length(vn1) != n_chain * n_draw) {
    koine_abort("netcdf_shape",
                "nchain * ndraw does not equal the number of draws (A1).")
  }
  vnames <- names(var_arrays)
  arr <- array(NA_real_, dim = c(n_draw, n_chain, length(vnames)),
               dimnames = list(NULL, NULL, vnames))
  for (k in seq_along(vnames)) {
    # var_arrays are (chain, draw); draws_array wants (draw, chain, var)
    arr[, , k] <- t(var_arrays[[vnames[k]]])
  }
  posterior::as_draws_array(arr)
}

# Reconstruct a provenance record from the attribute list.
.attrs_to_provenance <- function(attrs) {
  ik <- list(family = attrs$independence_family,
             kernel = attrs$independence_kernel)
  provenance_record(
    method = attrs$method, method_version = attrs$method_version,
    impl_language = attrs$impl_language, seed = attrs$seed,
    derived_seed = attrs$derived_seed, ir_hash = attrs$ir_hash,
    maturity_tier = attrs$maturity_tier, independence_key = ik,
    role = attrs$role, exactness = attrs$exactness,
    timestamp_iso8601 = attrs$timestamp_iso8601)
}
