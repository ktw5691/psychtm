#' @include model-class.R
NULL

#' An S4 class to represent a regression model.
#'
#' @slot ndocs The number of documents/observations.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot mu0 A (p + 1) x 1 matrix of prior means for eta.
#' @slot sigma0 A (p + 1) x (p + 1) prior covariance matrix for eta.
#' @slot a0 A prior shape hyperparameter for sigma2.
#' @slot b0 A prior rate hyperparameter for sigma2.
#' @slot eta_start A (p + 1) x 1 matrix of starting values for eta.
#' @slot eta A nchain x (p + 1) matrix of draws of regression coefficients.
#' @slot sigma2 A nchain x 1 numeric vector of draws of the residual variance.
#' @slot loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @slot logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot p_eff The effective number of parameters.
#' @slot lpd A nchain x ndocs matrix of predictive posterior likelihoods (NOT
#'   log-likelihoods).
Mlr <- setClass("Mlr",
  contains = "Model",
  slots = c(
    a0          = "numeric",
    b0          = "numeric",
    sigma2      = "numeric"),
  prototype = list(
    a0 = NA_real_,
    b0 = NA_real_,
    sigma2 = NA_real_
  )
)

#### sigma2
setGeneric("sigma2", function(x) standardGeneric("sigma2"))

setMethod("sigma2", "Mlr", function(x) x@sigma2)

setGeneric("sigma2<-", function(x, value) standardGeneric("sigma2<-"))

setMethod("sigma2<-", "Mlr", function(x, value) {
  x@sigma2 <- value
  x
})

#### a0
setGeneric("a0", function(x) standardGeneric("a0"))

setMethod("a0", "Mlr", function(x) x@a0)

setGeneric("a0<-", function(x, value) standardGeneric("a0<-"))

setMethod("a0<-", "Mlr", function(x, value) {
  x@a0 <- value
  x
})

#### b0
setGeneric("b0", function(x) standardGeneric("b0"))

setMethod("b0", "Mlr", function(x) x@b0)

setGeneric("b0<-", function(x, value) standardGeneric("b0<-"))

setMethod("b0<-", "Mlr", function(x, value) {
  x@b0 <- value
  x
})

# setMethod("initialize", "Mlr",
#           function(.Object) {
#             .Object <- callNextMethod(.Object)
#             .Object
#           }
# )
