#' An S4 class to represent a LDA model.
#'
#' @slot ntopics The number of topics for the LDA model (default: 2).
#' @slot ndocs The number of documents in the corpus.
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot alpha A numeric prior hyperparameter for theta.
#' @slot gamma A numeric prior hyperparameter for beta.
#' @slot loglike The log-likelihood (up to an additive constant).
#' @slot logpost The log-posterior (up to an additive constant).
Lda <- setClass("Lda",
                 slots = list(ntopics   = "numeric",
                              ndocs     = "numeric",
                              nvocab    = "numeric",
                              nchain    = "numeric",
                              topics    = "numeric",
                              alpha     = "numeric",
                              gamma     = "numeric",
                              loglike   = "numeric",
                              logpost   = "numeric"))

#' An S4 class to represent a sLDA model.
#'
#' @slot eta A M x K numeric matrix of draws of topic regression coefficients
#' @slot sigma2 A M x 1 numeric vector of draws of residual variance
#' @slot mu0 A K x 1 numeric matrix of prior means for eta
#' @slot sigma0 A K x K numeric prior covariance matrix for eta
#' @slot a0 A numeric prior shape hyperparameter for sigma2
#' @slot b0 A numeric prior rate hyperparameter for sigma2
#' @slot eta_start A K x 1 numeric matrix of starting values for eta
#' @slot ntopics The number of topics for the LDA model (default: 2).
#' @slot ndocs The number of documents in the corpus.
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot alpha A numeric prior hyperparameter for theta.
#' @slot gamma A numeric prior hyperparameter for beta.
#' @slot loglike The log-likelihood (up to an additive constant).
#' @slot logpost The log-posterior (up to an additive constant).
#' @slot p_eff The effective number of parameters.
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot lpd A matrix of iterations x observations of predictive posterior
#'   likelihoods (NOT log-likelihoods).
Slda <- setClass("Slda",
                 slots = list(eta       = "matrix",
                              sigma2    = "matrix",
                              mu0       = "numeric",
                              sigma0    = "matrix",
                              a0        = "numeric",
                              b0        = "numeric",
                              eta_start = "numeric",
                              p_eff     = "numeric",
                              waic      = "numeric",
                              se_waic   = "numeric",
                              lpd       = "matrix"),
                 contains = "Lda")

#' An S4 class to represent a sLDA logistic model.
#'
#' @slot eta A M x K numeric matrix of draws of topic regression coefficients
#' @slot mu0 A K x 1 numeric matrix of prior means for eta
#' @slot sigma0 A K x K numeric prior covariance matrix for eta
#' @slot eta_start A K x 1 numeric matrix of starting values for eta
#' @slot proposal_sd A K x 1 vector of proposal scales for Metropolis-Hastings
#'   sampling of eta.
#' @slot ntopics The number of topics for the LDA model (default: 2).
#' @slot ndocs The number of documents in the corpus.
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot alpha A numeric prior hyperparameter for theta.
#' @slot gamma A numeric prior hyperparameter for beta.
#' @slot loglike The log-likelihood (up to an additive constant).
#' @slot logpost The log-posterior (up to an additive constant).
#' @slot p_eff The effective number of parameters.
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot lpd A matrix of iterations x observations of predictive posterior
#'   likelihoods (NOT log-likelihoods).
Sldalogit <- setClass("Sldalogit",
                      slots = list(ndocs       = "numeric",
                                   nchain      = "numeric",
                                   eta         = "matrix",
                                   mu0         = "numeric",
                                   sigma0      = "matrix",
                                   eta_start   = "numeric",
                                   proposal_sd = "numeric",
                                   loglike     = "numeric",
                                   logpost     = "numeric",
                                   p_eff       = "numeric",
                                   waic        = "numeric",
                                   se_waic     = "numeric",
                                   lpd         = "matrix"),
                      contains = "Lda")

#' An S4 class to represent a logistic regression model.
#'
#' @slot eta A M x (p + 1) numeric matrix of draws of topic regression coefficients
#' @slot mu0 A (p + 1) x 1 numeric matrix of prior means for eta
#' @slot sigma0 A (p + 1) x (p + 1) numeric prior covariance matrix for eta
#' @slot eta_start A (p + 1) x 1 numeric matrix of starting values for eta
#' @slot proposal_sd A K x 1 vector of proposal scales for Metropolis-Hastings
#'   sampling of eta.
#' @slot ndocs The number of documents in the corpus.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot loglike The log-likelihood (up to an additive constant).
#' @slot logpost The log-posterior (up to an additive constant).
#' @slot p_eff The effective number of parameters.
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot lpd A matrix of iterations x observations of predictive posterior
#'   likelihoods (NOT log-likelihoods).
Logistic <- setClass("Logistic",
                      slots = list(nchain      = "numeric",
                                   ndocs       = "numeric",
                                   eta         = "matrix",
                                   mu0         = "numeric",
                                   sigma0      = "matrix",
                                   eta_start   = "numeric",
                                   proposal_sd = "numeric",
                                   loglike     = "numeric",
                                   logpost     = "numeric",
                                   p_eff       = "numeric",
                                   waic        = "numeric",
                                   se_waic     = "numeric",
                                   lpd         = "matrix"))

#' An S4 class to represent a regression model.
#'
#' @slot eta A M x (p + 1) numeric matrix of draws of topic regression
#'   coefficients
#' @slot sigma2 A M x 1 numeric vector of draws of residual variance
#' @slot mu0 A (p + 1) x 1 numeric matrix of prior means for eta
#' @slot sigma0 A (p + 1) x (p + 1) numeric prior covariance matrix for eta
#' @slot a0 A numeric prior shape hyperparameter for sigma2
#' @slot b0 A numeric prior rate hyperparameter for sigma2
#' @slot eta_start A (p + 1) x 1 numeric matrix of starting values for eta
#' @slot ndocs The number of documents in the corpus.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot loglike The log-likelihood (up to an additive constant).
#' @slot logpost The log-posterior (up to an additive constant).
#' @slot p_eff The effective number of parameters.
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot lpd A matrix of iterations x observations of predictive posterior
#'   likelihoods (NOT log-likelihoods).
Mlr <- setClass("Mlr",
                     slots = list(nchain      = "numeric",
                                  ndocs       = "numeric",
                                  eta         = "matrix",
                                  sigma2      = "numeric",
                                  mu0         = "numeric",
                                  sigma0      = "matrix",
                                  a0          = "numeric",
                                  b0          = "numeric",
                                  eta_start   = "numeric",
                                  loglike     = "numeric",
                                  logpost     = "numeric",
                                  p_eff       = "numeric",
                                  waic        = "numeric",
                                  se_waic     = "numeric",
                                  lpd         = "matrix"))

setMethod("initialize", "Lda",
          function(.Object, alpha = 0.1, gamma = 1.01, a0 = 0.001, b0 = 0.001) {
            .Object <- callNextMethod(.Object)
            .Object@ntopics <- 2
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

setMethod("initialize", "Lda",
          function(.Object, alpha = 0.1, gamma = 1.01) {
            .Object <- callNextMethod(.Object)
            .Object@ntopics <- 2
            .Object@alpha <- alpha
            .Object@gamma <- gamma
            .Object@loglike = NaN
            .Object@logpost = NaN
            .Object
          }
)

setMethod("initialize", "Logistic",
          function(.Object, proposal_sd = 2.38) {
            .Object <- callNextMethod(.Object)
            .Object
          }
)


setMethod("initialize", "Sldalogit",
          function(.Object, alpha = 0.1, gamma = 1.01, proposal_sd = 2.38) {
            .Object <- callNextMethod(.Object)
            .Object
          }
)
