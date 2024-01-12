

#' Internal function doing standard checks
#'
#' @param St3_observed numeric: 0-1
#' @param St2_observed numeric: 0-1
#' @param St1_observed numeric: 0-1
#' @param St3_true numeric: 0-1
#' @param St2_true numeric: 0-1
#' @param St1_true numeric: 0-1
#' @param theta_01 numeric
#' @param theta_02 numeric
#' @param err numeric
#' @param mu numeric
#'
#' @return error if necessary
#' @keywords internal
sme_checks <- function(
  St3_observed,
  St2_observed,
  St1_observed,
  St3_true,
  St2_true,
  St1_true,
  theta_01,
  theta_02,
  err,
  mu,
  only_2 = TRUE
){

  if (only_2 == TRUE) {
    if (!(St3_observed == 0 || St3_observed == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }
    if (!(St3_true == 0     || St3_true == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }
  }
    if (!(St2_observed == 0 || St2_observed == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }
    if (!(St1_observed == 0 || St1_observed == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }
    if (!(St2_true == 0     || St2_true == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }
    if (!(St1_true == 0     || St1_true == 1)) {
      stop("St3_observed contains values not 0 or 1")
    }

  if (!is.numeric(theta_01)) stop("argument `theta_01` must be numeric")
  if (!is.numeric(theta_02)) stop("argument `theta_02` must be numeric")
  if (!is.numeric(err))      stop("argument `err` must be numeric")
  if (!is.numeric(mu))       stop("argument `mu` must be numeric")
}








