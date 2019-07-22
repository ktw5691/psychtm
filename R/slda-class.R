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
                slots = list(ndocs       = "numeric",
                             nchain      = "numeric",
                             mu0         = "numeric",
                             sigma0      = "matrix",
                             a0          = "numeric",
                             b0          = "numeric",
                             eta_start   = "numeric",
                             eta         = "matrix",
                             sigma2      = "numeric",
                             loglike     = "numeric",
                             logpost     = "numeric",
                             waic        = "numeric",
                             se_waic     = "numeric",
                             p_eff       = "numeric",
                             lpd         = "matrix"))

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
                      slots = list(ndocs       = "numeric",
                                   nchain      = "numeric",
                                   mu0         = "numeric",
                                   sigma0      = "matrix",
                                   eta_start   = "numeric",
                                   eta         = "matrix",
                                   proposal_sd = "numeric",
                                   loglike     = "numeric",
                                   logpost     = "numeric",
                                   waic        = "numeric",
                                   se_waic     = "numeric",
                                   p_eff       = "numeric",
                                   lpd         = "matrix"))

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
                  slots = list(ndocs       = "numeric",
                               nvocab      = "numeric",
                               ntopics     = "numeric",
                               nchain      = "numeric",
                               alpha       = "numeric",
                               gamma       = "numeric",
                               mu0         = "numeric",
                               sigma0      = "matrix",
                               a0          = "numeric",
                               b0          = "numeric",
                               eta_start   = "numeric",
                               topics      = "array",
                               eta         = "matrix",
                               sigma2      = "numeric",
                               proposal_sd = "numeric",
                               loglike     = "numeric",
                               logpost     = "numeric",
                               waic        = "numeric",
                               se_waic     = "numeric",
                               p_eff       = "numeric",
                               lpd         = "matrix"))

setMethod("initialize", "Mlr",
          function(.Object) {
            .Object <- callNextMethod(.Object)
            .Object
          }
)

setMethod("initialize", "Logistic",
          function(.Object) {
            .Object <- callNextMethod(.Object)
            .Object
          }
)

setMethod("initialize", "Sldax",
          function(.Object, alpha = 0.1, gamma = 1.01, a0 = 0.001, b0 = 0.001) {
            .Object <- callNextMethod(.Object)
            .Object@ntopics <- 2
            .Object@alpha   <- alpha
            .Object@gamma   <- gamma
            .Object@a0      <- a0
            .Object@b0      <- b0
            .Object@loglike <- NaN
            .Object@logpost <- NaN
            .Object
          }
)
