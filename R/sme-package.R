#' sme: Functions to use structural misclassification estimator
#'
#' This package allows the use of different estimators to identify
#' dynamics of discrete variables that are subject to classification error
#'
#' @section sme functions: The sme functions have certain distinctions
#'   First, `ar1` refers to the fact that the true binary process is an AR1
#'   process possessing the Markov property, while `ar2` is is truly AR2 and
#'   does not possess the Markov property.
#'   Second, `sym` or `asym` refers to whether classification error is symmetric
#'   around the true status or not.
#'   Third, `tenure` is highly specialized and should be used only if you
#'   are very confident of your specification. It refers to the use of a
#'   tenure/duration variable that can be included in the structural equations
#'   to identify classification error.
#'   Fourth, `inc` refers to a more ad-hoc use of errors in covariates than
#'   found in `tenure`, where it exploits exogenously determined inconsistencies
#'   in observed covariates to identify classification error
#'
#' @docType package
#' @name sme
#' @import collapse
#' @import mlrMBO
#' @import maxLik

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE

# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c(
      ".",
      ".I",
      ".N",
      ".SD",
      ".",
      "!!",
      ":="
    ),
    package = utils::packageName()
  )
}


NULL
