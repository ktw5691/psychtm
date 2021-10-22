#' psychtm: A package for text mining methods for psychological research
#'
#' The `psychtm` package provides estimation, summarization, and goodness-of-fit
#' functions:
#'
#' @section Model Fitting:
#' The workhorse function for Bayesian estimation of topic models is
#' [`gibbs_sldax()`][gibbs_sldax()].
#' Similarly, see [`gibbs_mlr()`][gibbs_mlr()] and
#' [`gibbs_logistic()`][gibbs_logistic()] to estimate regression models with
#' continuous and dichotomous outcomes, respectively.
#'
#' @section Parameter Estimates and Goodness-of-Fit:
#' See [`sldax-summary`] for functions to obtain and summarize parameter
#' estimates and to compute goodness-of-fit metrics.
#'
#' @docType package
#' @name psychtm
NULL
#> NULL

#' @importFrom methods setClass setGeneric setMethod setRefClass is new validObject
NULL

#' @importFrom stats median model.matrix model.response quantile
NULL
