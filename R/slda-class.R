#' An S4 class to represent a LDA model.
#'
#' @slot beta A K x V x M numeric array of draws of topic-word probabilities
#' @slot theta A D x K x M numeric array of draws of document-topic
#' probabilities
#' @slot alpha A numeric prior hyperparameter for theta
#' @slot gamma A numeric prior hyperparameter for beta
Lda <- setClass("Lda",
                 slots = list(ntopics   = "numeric",
                              ndocs     = "numeric",
                              nvocab    = "numeric",
                              nchain    = "numeric",
                              beta      = "array",
                              theta     = "array",
                              alpha     = "numeric",
                              gamma     = "numeric",
                              loglike   = "numeric",
                              logpost   = "numeric"))

#' An S4 class to represent a sLDA model.
#'
#' @slot eta A M x K numeric matrix of draws of topic regression coefficients
#' @slot sigma2 A M x 1 numeric vector of draws of residual variance
#' @slot beta A K x V x M numeric array of draws of topic-word probabilities
#' @slot theta A D x K x M numeric array of draws of document-topic
#' probabilities
#' @slot mu0 A K x 1 numeric matrix of prior means for eta
#' @slot sigma0 A K x K numeric prior covariance matrix for eta
#' @slot alpha A numeric prior hyperparameter for theta
#' @slot gamma A numeric prior hyperparameter for beta
#' @slot a0 A numeric prior shape hyperparameter for sigma2
#' @slot b0 A numeric prior rate hyperparameter for sigma2
#' @slot eta_start A K x 1 numeric matrix of starting values for eta
Slda <- setClass("Slda",
  slots = list(eta       = "matrix",
               sigma2    = "matrix",
               mu0       = "matrix",
               sigma0    = "matrix",
               a0        = "numeric",
               b0        = "numeric",
               eta_start = "matrix"),
  contains = "Lda")

setMethod("initialize", "Lda",
          function(.Object, alpha = 0.1, gamma = 1.01, a0 = 0.001, b0 = 0.001) {
            .Object <- callNextMethod(.Object)
            .Object@alpha <- alpha
            .Object@gamma <- gamma
            .Object@loglike = NaN
            .Object@logpost = NaN
            .Object
          }
)

setMethod("initialize", "Slda",
          function(.Object, alpha = 0.1, gamma = 1.01, a0 = 0.001, b0 = 0.001) {
            .Object <- callNextMethod(.Object)
            .Object@a0 = a0
            .Object@b0 = b0
            .Object
          }
)
