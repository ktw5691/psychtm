#' An S4 super class to represent a regression-like model
#'
#' @slot ndocs The number of documents/observations.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot mu0 A (p + 1) x 1 matrix of prior means for eta.
#' @slot sigma0 A (p + 1) x (p + 1) prior covariance matrix for eta.
#' @slot eta_start A (p + 1) x 1 matrix of starting values for eta.
#' @slot eta A nchain x (p + 1) matrix of draws of regression coefficients.
#' @slot loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @slot logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot p_eff The effective number of parameters.
#' @slot lpd A nchain x ndocs matrix of predictive posterior likelihoods.
#' @slot extra A list of additional model fitting information. Contains
#'   time_elapsed, start_time, end_time, corrected_label_switching, and call.
#'
#' @param ndocs The number of documents/observations.
#' @param nchain The number of iterations of the Gibbs sampler.
#' @param mu0 A (p + 1) x 1 matrix of prior means for eta.
#' @param sigma0 A (p + 1) x (p + 1) prior covariance matrix for eta.
#' @param eta_start A (p + 1) x 1 matrix of starting values for eta.
#' @param eta A nchain x (p + 1) matrix of draws of regression coefficients.
#' @param loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @param logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @param waic WAIC (up to an additive constant) on the deviance scale.
#' @param se_waic Standard error of the WAIC.
#' @param p_eff The effective number of parameters.
#' @param lpd A nchain x ndocs matrix of predictive posterior likelihoods.
#' @param x An `Model` object.
#' @param value A value to assign to a slot for `x`
#'
#' @name Model-class
#' @rdname Model-class
#' @keywords classes
#' @exportClass Model
setClass("Model",
         slots = c(
           ndocs       = "numeric",
           nchain      = "numeric",
           mu0         = "numeric",
           sigma0      = "matrix",
           eta_start   = "numeric",
           eta         = "matrix",
           loglike     = "numeric",
           logpost     = "numeric",
           waic        = "numeric",
           se_waic     = "numeric",
           p_eff       = "numeric",
           lpd         = "matrix",
           extra       = "list"),
         prototype = list(
           ndocs = NA_real_,
           nchain = NA_real_,
           mu0 = NA_real_,
           sigma0 = matrix(NA_real_),
           eta_start = NA_real_,
           eta = matrix(NA_real_),
           loglike = NA_real_,
           logpost = NA_real_,
           waic = NA_real_,
           se_waic = NA_real_,
           p_eff = NA_real_,
           lpd = matrix(NA_real_),
           extra = list(time_elapsed = NA_real_,
                        start_time = NA_real_,
                        end_time = NA_real_,
                        corrected_label_switching = FALSE,
                        call = NA_character_)
         )
)

#' S4 class for a regression model that inherits from [Model-class].
#'
#' @slot a0 A prior shape hyperparameter for sigma2.
#' @slot b0 A prior rate hyperparameter for sigma2.
#' @slot sigma2 A nchain x 1 numeric vector of draws of the residual variance.
#'
#' @param a0 A prior shape hyperparameter for sigma2.
#' @param b0 A prior rate hyperparameter for sigma2.
#' @param sigma2 A nchain x 1 numeric vector of draws of the residual variance.
#' @param x An `Model` object.
#' @param value A value to assign to a slot for `x`
#' @param ...	additional arguments to be passed to the low level regression
#'   fitting functions (see below).
#'
#' @name Mlr-class
#' @rdname Mlr-class
#' @keywords classes
#' @exportClass Mlr
setClass("Mlr",
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

#' S4 class for a logistic regression model that inherits from [Model-class]
#'
#' @slot proposal_sd A vector of p + 1 proposal scales/standard deviations for
#'   sampling of p + 1 regression coefficients by Metropolis-Hastings.
#' @param proposal_sd A vector of p + 1 proposal scales/standard deviations for
#'   sampling of p + 1 regression coefficients by Metropolis-Hastings.
#' @param x An `Logistic` object.
#' @param value A value to assign to a slot for `x`
#' @param ...	additional arguments to be passed to the low level regression
#'   fitting functions (see below).
#'
#' @name Logistic-class
#' @rdname Logistic-class
#' @keywords classes
#' @exportClass Logistic
setClass("Logistic",
         contains = "Model",
         slots = c(
           proposal_sd = "numeric"
         ),
         prototype = list(
           proposal_sd = NA_real_
         )
)

#' S4 class to represent a SLDAX general model that inherits from [Mlr-class]
#'   and [Logistic-class].
#'
#' @slot nvocab The number of terms in the corpus vocabulary.
#' @param nvocab The number of terms in the corpus vocabulary.
#' @slot ntopics The number of topics for the LDA model.
#' @param ntopics The number of topics for the LDA model (default: `2`).
#' @slot alpha A numeric prior hyperparameter for theta.
#' @param alpha A numeric prior hyperparameter for theta (default: `1.0`).
#' @slot gamma A numeric prior hyperparameter for beta.
#' @param gamma A numeric prior hyperparameter for beta  (default: `1.0`).
#' @slot topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @param topics A D x max(N_d) x M numeric array of topic draws. 0 indicates an
#'   unused word index (i.e., the document did not have a word at that index).
#' @slot theta A D x K x M numeric array of topic proportions.
#' @param theta A D x K x M numeric array of topic proportions.
#' @slot beta A K x V x M numeric array of topic-vocabulary distributions.
#' @param beta A K x V x M numeric array of topic-vocabulary distributions.
#' @param x An `Sldax` object.
#' @param value A value to assign to a slot for `x`
#' @param ...	additional arguments to be passed to the low level regression
#'   fitting functions (see below).
#'
#' @name Sldax-class
#' @rdname Sldax-class
#' @keywords classes
#' @exportClass Sldax
setClass("Sldax",
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
