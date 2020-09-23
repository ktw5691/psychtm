#' @include mlr-class.R logistic-class.R
NULL

#' An S4 class to represent a sLDAX general model.
#'
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @slot ntopics The number of topics for the LDA model (default: 2).
#' @slot alpha A numeric prior hyperparameter for theta.
#' @slot gamma A numeric prior hyperparameter for beta.
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot theta A D x K x M numeric array of topic proportions.
#' @slot beta A K x V x M numeric array of topic-vocabulary distributions.
Sldax <- setClass("Sldax",
  contains = c("Mlr", "Logistic"),
  slots = c(
    nvocab      = "numeric",
    ntopics     = "numeric",
    alpha       = "numeric",
    gamma       = "numeric",
    topics      = "array",
    theta       = "array",
    beta        = "array"),
  prototype = list(
    nvocab = NA_real_,
    ntopics = NA_real_,
    alpha = NA_real_,
    gamma = NA_real_,
    topics = array(NA_real_, dim = c(1, 2, 1)),
    theta = array(NA_real_, dim = c(1, 2, 1)),
    beta = array(NA_real_, dim = c(2, 2, 1))
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

#### theta
setGeneric("theta", function(x) standardGeneric("theta"))

setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))

setMethod("theta", "Sldax", function(x) x@theta)

setMethod("theta<-", "Sldax", function(x, value) {
  x@theta <- value
  x
})

#### beta
setGeneric("beta_", function(x) standardGeneric("beta_"))

setGeneric("beta_<-", function(x, value) standardGeneric("beta_<-"))

setMethod("beta_", "Sldax", function(x) x@beta)

setMethod("beta_<-", "Sldax", function(x, value) {
  x@beta <- value
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

#' Helper function (constructor) for Sldax class.
#'
Sldax <- function(nvocab, topics, theta, beta, ntopics = 2.0, alpha = 0.1,
                  gamma = 1.01, ...) {
  super <- Model(...)
  nvocab <- as.double(nvocab)
  topics <- as.array(topics)
  theta <- as.array(theta)
  beta <- as.array(beta)
  ntopics <- as.double(ntopics)
  alpha <- as.double(alpha)
  gamma <- as.double(gamma)

  new("Sldax", nvocab = nvocab, topics = topics, theta = theta, beta = beta,
      ntopics = ntopics, alpha = alpha, gamma = gamma,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

#' Validator function for Sldax class
#' TODO
