#' @keywords internal
#' @aliases koine-package
"_PACKAGE"

## ----------------------------------------------------------------------------
## Package-level notes
##
## koine synthesises one independent predictive opinion for a Bayesian model
## from a battery of inference lenses, or abstains. The contract surface this
## package implements is specified in `koine-stage0-contracts-v1.md` (Stage-0).
## Section references in the roxygen blocks (e.g. "SS5", "F1", "A3") point at
## that document.
##
## Two cardinal invariants hold everywhere (Stage-0 SS1):
##   I1  Honest representation  --- never fabricate a value not computed;
##                                  absence is an explicit NULL plus a reason.
##   I2  Predictive space only  --- combination happens on p(y-tilde | .); never
##                                  average parameter posteriors across
##                                  heterogeneous structures.
## ----------------------------------------------------------------------------

# The NIMBLE lens writes its model in NIMBLE's BUGS dialect inside nimbleCode();
# the tokens below are NIMBLE model-language symbols, not R objects, so they are
# declared here to keep `R CMD check` from flagging them as undefined globals.
utils::globalVariables(c("K", "N", "X", "log_sigma", "inprod",
                         "log<-", "logit<-"))
