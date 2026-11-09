#' @keywords internal
"_PACKAGE"

# The indicator functions use data.table's `:=` and `[.data.table` with
# non-standard evaluation. data.table only allows that from packages it
# considers "data.table-aware" (see vignette("datatable-importing")); without
# this flag those calls error with "[ was called on a data.table in an
# environment that is not data.table-aware".
.datatable.aware <- TRUE

#' @importFrom data.table :=
#' @importFrom utils globalVariables
NULL

# data.table NSE symbols and column names used inside `[.data.table` calls
utils::globalVariables(c(".N", ".SD", "abund", "area", "cum_abund", "cum_area"))
