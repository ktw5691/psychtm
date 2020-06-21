#' @include mlr-class.R logistic-class.R
NULL

#' An S4 class to represent a sLDAX general model.
#'
#' @slot ndocs The number of documents in the corpus.
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @slot ntopics The number of topics for the LDA model (default: 2).
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot alpha A numeric prior hyperparameter for theta.
#' @slot gamma A numeric prior hyperparameter for beta.
#' @slot mu0 A q x 1 numeric matrix of prior means for eta.
#' @slot sigma0 A q x q numeric prior covariance matrix for eta.
#' @slot a0 A numeric prior shape hyperparameter for sigma2.
#' @slot b0 A numeric prior rate hyperparameter for sigma2.
#' @slot eta_start A q x 1 numeric matrix of starting values for eta.
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot eta A nchain x q numeric matrix of draws of topic regression coefficients.
#' @slot sigma2 A nchain x 1 numeric vector of draws of residual variance.
#' @slot proposal_sd A q x 1 vector of proposal scales for Metropolis-Hastings
#'   sampling of eta.
#' @slot loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @slot logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot p_eff The effective number of parameters.
#' @slot lpd A nchain x ndocs matrix of predictive posterior likelihoods (NOT
#'   log-likelihoods).
Sldax <- setClass("Sldax",
  contains = c("Mlr", "Logistic"),
  slots = c(
    nvocab      = "numeric",
    ntopics     = "numeric",
    alpha       = "numeric",
    gamma       = "numeric",
    topics      = "array"),
  prototype = list(
    nvocab = NA_real_,
    ntopics = NA_real_,
    alpha = NA_real_,
    gamma = NA_real_,
    topics = array(NA_real_, dim = c(1, 1, 1))
  )
)

#### topics
setGeneric("topics", function(x) standardGeneric("topics"))

setGeneric("topics<-", function(x, value) standardGeneric("topics<-"))

setMethod("topics", "Sldax", function(x) x@topics)

setMethod("topics<-", "Sldax", function(x, value) {
  x@topics <- value
  x
})

#### gamma
setGeneric("gamma_", function(x) standardGeneric("gamma_"))
setGeneric("gamma_<-", function(x, value) standardGeneric("gamma_<-"))
setMethod("gamma_", "Sldax", function(x) x@gamma)
setMethod("gamma_<-", "Sldax", function(x, value) {
  x@gamma <- value
  x
})

#### alpha
setGeneric("alpha", function(x) standardGeneric("alpha"))
setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))
setMethod("alpha", "Sldax", function(x) x@alpha)
setMethod("alpha<-", "Sldax", function(x, value) {
  x@alpha <- value
  x
})

#### ntopics
setGeneric("ntopics", function(x) standardGeneric("ntopics"))

setGeneric("ntopics<-", function(x, value) standardGeneric("ntopics<-"))

setMethod("ntopics", "Sldax", function(x) x@ntopics)

setMethod("ntopics<-", "Sldax", function(x, value) {
  x@ntopics <- value
  x
})

#### nvocab
setGeneric("nvocab", function(x) standardGeneric("nvocab"))

setGeneric("nvocab<-", function(x, value) standardGeneric("nvocab<-"))

setMethod("nvocab", "Sldax", function(x) x@nvocab)

setMethod("nvocab<-", "Sldax", function(x, value) {
  x@nvocab <- value
  x
})
#
# setMethod("initialize", "Sldax",
#           function(.Object, alpha = 0.1, gamma = 1.01, a0 = 0.001, b0 = 0.001) {
#             .Object <- callNextMethod(.Object)
#             .Object@ntopics <- 2
#             .Object@alpha   <- alpha
#             .Object@gamma   <- gamma
#             .Object@a0      <- a0
#             .Object@b0      <- b0
#             .Object@loglike <- NaN
#             .Object@logpost <- NaN
#             .Object
#           }
# )
