#' @keywords internal
"_PACKAGE"

#' @importFrom methods is
#' @importFrom stats median p.adjust sd setNames ks.test fisher.test
#' @importFrom utils URLencode data stack
NULL

# Data sets are loaded into the calling frame via utils::data(), which the
# code analyser cannot see. Declare them so R CMD check does not report them
# as undefined globals.
utils::globalVariables(c("NCBI_ids", "TaxSEA_db"))
