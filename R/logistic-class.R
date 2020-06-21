#' @include model-class.R
NULL

#' An S4 class to represent a logistic regression model.
#'
#' @slot ndocs The number of documents/observations.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot mu0 A (p + 1) x 1 matrix of prior means for eta.
#' @slot sigma0 A (p + 1) x (p + 1) prior covariance matrix for eta.
#' @slot eta_start A (p + 1) x 1 matrix of starting values for eta.
#' @slot eta A nchain x (p + 1) matrix of draws of regression coefficients.
#' @slot proposal_sd A (p + 1) x 1 vector of proposal scales for
#'   Metropolis-Hastings sampling of eta.
#' @slot loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @slot logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot p_eff The effective number of parameters.
#' @slot lpd A nchain x ndocs matrix of predictive posterior likelihoods (NOT
#'   log-likelihoods).
Logistic <- setClass("Logistic",
  contains = "Model",
  slots = c(
    proposal_sd = "numeric"
  ),
  prototype = list(
    proposal_sd = NA_real_
  )
)

#### proposal_sd
setGeneric("proposal_sd", function(x) standardGeneric("proposal_sd"))

setGeneric("proposal_sd<-", function(x, value) standardGeneric("proposal_sd<-"))

setMethod("proposal_sd", "Logistic", function(x) x@proposal_sd)

setMethod("proposal_sd<-", "Logistic", function(x, value) {
  x@proposal_sd <- value
  x
})

# setMethod("initialize", "Logistic",
#           function(.Object) {
#             .Object <- callNextMethod(.Object)
#             .Object
#           }
# )
