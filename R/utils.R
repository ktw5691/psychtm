#' Check for integer argument
#'
#' @param x Argument to check.
#' @param tol Tolerance (default: `sqrt(.Machine$double.eps)`)
#'
#' @noRd
is.whole_number <- function(x, tol = sqrt(.Machine$double.eps)) {
  return(abs(x - round(x)) < tol)
}

#' Check for non-negative integer
#'
#' @param x Argument to check.
#' @param tol Tolerance (default: `sqrt(.Machine$double.eps)`)
#'
#' @noRd
is.non_negative_integer <- function(x, tol = sqrt(.Machine$double.eps)) {
  return(abs(x - round(x)) < tol & x >= 0)
}

#' Check for positive integer
#'
#' @param x Argument to check.
#' @param tol Tolerance (default: `sqrt(.Machine$double.eps)`)
#'
#' @noRd
is.positive_integer <- function(x, tol = sqrt(.Machine$double.eps)) {
  return(abs(x - round(x)) < tol & x > 0)
}

#' Check for missing argument
#'
#' @param x Argument to check.
#'
#' @noRd
missing_msg <- function(x) {
  stop(paste0('\"', x, '\" argument is missing.'), call. = FALSE)
}

#' Check for logical and non-missing argument
#'
#' @param arg Argument to check.
#'
#' @noRd
check_logical <- function(arg) {
  good <- FALSE
  if (is.logical(arg) & !is.na(arg)) good <- TRUE
  return(good)
}
