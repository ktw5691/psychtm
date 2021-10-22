########----------------------------------------------------------------########
########                            Model Class                         ########
########----------------------------------------------------------------########

# ndocs
#' @rdname Model
#' @export
setGeneric("ndocs", function(x) standardGeneric("ndocs"))

#' @rdname Model
setMethod("ndocs", "Model", function(x) x@ndocs)

#' @rdname Model
#' @export
setGeneric("ndocs<-", function(x, value) standardGeneric("ndocs<-"))

#' @rdname Model
setMethod("ndocs<-", "Model", function(x, value) {
  x@ndocs <- value
  x
})

# nchain
#' @rdname Model
#' @export
setGeneric("nchain", function(x) standardGeneric("nchain"))

#' @rdname Model
setMethod("nchain", "Model", function(x) x@nchain)

#' @rdname Model
#' @export
setGeneric("nchain<-", function(x, value) standardGeneric("nchain<-"))

#' @rdname Model
setMethod("nchain<-", "Model", function(x, value) {
  x@nchain <- value
  x
})

# mu0
#' @rdname Model
#' @export
setGeneric("mu0", function(x) standardGeneric("mu0"))

#' @rdname Model
setMethod("mu0", "Model", function(x) x@mu0)

#' @rdname Model
#' @export
setGeneric("mu0<-", function(x, value) standardGeneric("mu0<-"))

#' @rdname Model
setMethod("mu0<-", "Model", function(x, value) {
  x@mu0 <- value
  x
})

# sigma0
#' @rdname Model
#' @export
setGeneric("sigma0", function(x) standardGeneric("sigma0"))

#' @rdname Model
setMethod("sigma0", "Model", function(x) x@sigma0)

#' @rdname Model
#' @export
setGeneric("sigma0<-", function(x, value) standardGeneric("sigma0<-"))

#' @rdname Model
setMethod("sigma0<-", "Model", function(x, value) {
  x@sigma0 <- value
  x
})

# eta_start
#' @rdname Model
#' @export
setGeneric("eta_start", function(x) standardGeneric("eta_start"))

#' @rdname Model
setMethod("eta_start", "Model", function(x) x@eta_start)

#' @rdname Model
#' @export
setGeneric("eta_start<-", function(x, value) standardGeneric("eta_start<-"))

#' @rdname Model
setMethod("eta_start<-", "Model", function(x, value) {
  x@eta_start <- value
  x
})

# eta
#' @rdname Model
#' @export
setGeneric("eta", function(x) standardGeneric("eta"))

#' @rdname Model
setMethod("eta", "Model", function(x) x@eta)

#' @rdname Model
#' @export
setGeneric("eta<-", function(x, value) standardGeneric("eta<-"))

#' @rdname Model
setMethod("eta<-", "Model", function(x, value) {
  x@eta <- value
  x
})

# loglike
#' @rdname Model
#' @export
setGeneric("loglike", function(x) standardGeneric("loglike"))

#' @rdname Model
setMethod("loglike", "Model", function(x) x@loglike)

#' @rdname Model
#' @export
setGeneric("loglike<-", function(x, value) standardGeneric("loglike<-"))

#' @rdname Model
setMethod("loglike<-", "Model", function(x, value) {
  x@loglike <- value
  x
})

# logpost
#' @rdname Model
#' @export
setGeneric("logpost", function(x) standardGeneric("logpost"))

#' @rdname Model
setMethod("logpost", "Model", function(x) x@logpost)

#' @rdname Model
#' @export
setGeneric("logpost<-", function(x, value) standardGeneric("logpost<-"))

#' @rdname Model
setMethod("logpost<-", "Model", function(x, value) {
  x@logpost <- value
  x
})

# waic
#' @rdname Model
#' @export
setGeneric("waic", function(x) standardGeneric("waic"))

#' @rdname Model
setMethod("waic", "Model", function(x) x@waic)

#' @rdname Model
#' @export
setGeneric("waic<-", function(x, value) standardGeneric("waic<-"))

#' @rdname Model
setMethod("waic<-", "Model", function(x, value) {
  x@waic <- value
  x
})

# se_waic
#' @rdname Model
#' @export
setGeneric("se_waic", function(x) standardGeneric("se_waic"))

#' @rdname Model
setMethod("se_waic", "Model", function(x) x@se_waic)

#' @rdname Model
#' @export
setGeneric("se_waic<-", function(x, value) standardGeneric("se_waic<-"))

#' @rdname Model
setMethod("se_waic<-", "Model", function(x, value) {
  x@se_waic <- value
  x
})

# p_eff
#' @rdname Model
#' @export
setGeneric("p_eff", function(x) standardGeneric("p_eff"))

#' @rdname Model
setMethod("p_eff", "Model", function(x) x@p_eff)

#' @rdname Model
#' @export
setGeneric("p_eff<-", function(x, value) standardGeneric("p_eff<-"))

#' @rdname Model
setMethod("p_eff<-", "Model", function(x, value) {
  x@p_eff <- value
  x
})

# lpd
#' @rdname Model
#' @export
setGeneric("lpd", function(x) standardGeneric("lpd"))

#' @rdname Model
setMethod("lpd", "Model", function(x) x@lpd)

#' @rdname Model
#' @export
setGeneric("lpd<-", function(x, value) standardGeneric("lpd<-"))

#' @rdname Model
setMethod("lpd<-", "Model", function(x, value) {
  x@lpd <- value
  x
})

# extra
#' @rdname Model
#' @export
setGeneric("extra", function(x) standardGeneric("extra"))

#' @rdname Model
setMethod("extra", "Model", function(x) x@extra)

#' @rdname Model
#' @export
setGeneric("extra<-", function(x, value) standardGeneric("extra<-"))

#' @rdname Model
setMethod("extra<-", "Model", function(x, value) {
  x@extra <- value
  x
})

# Helper function (constructor) for Model class.
#' @rdname Model
Model <- function(ndocs, nchain = 1, mu0 = NaN, sigma0 = NaN,
                  eta_start = NaN, eta = NaN, loglike = NaN, logpost = NaN,
                  waic = NaN, se_waic = NaN, p_eff = NaN, lpd = NaN) {
  ndocs <- as.double(ndocs)
  nchain <- as.double(nchain)
  mu0 <- as.double(mu0)
  sigma0 <- as.matrix(sigma0)
  eta_start <- as.double(mu0)
  eta <- as.matrix(eta)
  loglike <- as.double(loglike)
  logpost <- as.double(logpost)
  waic <- as.double(waic)
  se_waic <- as.double(se_waic)
  p_eff <- as.double(p_eff)
  lpd <- as.matrix(lpd)

  new("Model", ndocs = ndocs, nchain = nchain, mu0 = mu0, sigma0 = sigma0,
      eta_start = eta_start, eta = eta, loglike = loglike,
      logpost = logpost, waic = waic, se_waic = se_waic, p_eff = p_eff,
      lpd = lpd)
}

# Validator function for Model class
setValidity("Model", function(object) {
  if ( (object@nchain != NROW(object@eta)) |
       (object@nchain != length(object@loglike)) |
       (object@nchain != length(object@logpost)) |
       (object@nchain != NROW(object@lpd))) {
    "@eta and @lpd must have the same number of rows as the length of @nchain"
  } else if ( (length(object@waic) != length(object@se_waic)) |
              (length(object@waic) != length(object@p_eff))) {
    "@waic, @se_waic, and @p_eff should all be of length 1"
  } else if ( (ncol(object@lpd) != object@ndocs)) {
    "@lpd should have number of columns equal to value of @ndocs"
  } else if ( (length(object@mu0) != length(object@eta_start)) |
              (length(object@mu0) != ncol(object@eta)) |
              (length(object@mu0) != nrow(object@sigma0)) |
              (length(object@mu0) != ncol(object@sigma0))) {
    "@mu0, @eta_start, and @eta should all have the same length which should be equal to the number of rows and columns in @sigma0"
  } else {
    TRUE
  }
})

########----------------------------------------------------------------########
########                             Mlr Class                          ########
########----------------------------------------------------------------########

#' @rdname Mlr
setGeneric("sigma2", function(x) standardGeneric("sigma2"))

#' @rdname Mlr
setMethod("sigma2", "Mlr", function(x) x@sigma2)

#' @rdname Mlr
setGeneric("sigma2<-", function(x, value) standardGeneric("sigma2<-"))

#' @rdname Mlr
setMethod("sigma2<-", "Mlr", function(x, value) {
  x@sigma2 <- value
  x
})

#' @rdname Mlr
setGeneric("a0", function(x) standardGeneric("a0"))

#' @rdname Mlr
setMethod("a0", "Mlr", function(x) x@a0)

#' @rdname Mlr
setGeneric("a0<-", function(x, value) standardGeneric("a0<-"))

#' @rdname Mlr
setMethod("a0<-", "Mlr", function(x, value) {
  x@a0 <- value
  x
})

#' @rdname Mlr
setGeneric("b0", function(x) standardGeneric("b0"))

#' @rdname Mlr
setMethod("b0", "Mlr", function(x) x@b0)

#' @rdname Mlr
setGeneric("b0<-", function(x, value) standardGeneric("b0<-"))

#' @rdname Mlr
setMethod("b0<-", "Mlr", function(x, value) {
  x@b0 <- value
  x
})

# Helper function (constructor) for Mlr class.
#' @rdname Mlr
Mlr <- function(a0 = 0.001, b0 = 0.001, sigma2 = NaN, ...) {
  super <- Model(...)
  a0 <- as.double(a0)
  b0 <- as.double(b0)
  sigma2 <- as.double(sigma2)

  new("Mlr", a0 = a0, b0 = b0, sigma2 = sigma2,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

########----------------------------------------------------------------########
########                          Logistic Class                        ########
########----------------------------------------------------------------########

#' Slot `@proposal_sd` generic
#' @rdname Logistic
setGeneric("proposal_sd", function(x) standardGeneric("proposal_sd"))

#' Slot `@proposal_sd` generic accessor
#' @rdname Logistic
setMethod("proposal_sd", "Logistic", function(x) x@proposal_sd)

#' Slot `@proposal_sd` generic setter function
#' @rdname Logistic
setGeneric("proposal_sd<-", function(x, value) standardGeneric("proposal_sd<-"))

#' Slot `@proposal_sd` setter method
#' @rdname Logistic
setMethod("proposal_sd<-", "Logistic", function(x, value) {
  x@proposal_sd <- value
  validObject(x)
  x
})

#' Logistic helper function (constructor) for `Logistic` class.
#' @rdname Logistic
Logistic <- function(proposal_sd = NaN, ...) {
  super <- Model(...)
  proposal_sd <- as.double(proposal_sd)

  new("Logistic", proposal_sd = proposal_sd,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

########----------------------------------------------------------------########
########                            Sldax Class                         ########
########----------------------------------------------------------------########

# topics
setGeneric("topics", function(x) standardGeneric("topics"))

#' @rdname Sldax
setMethod("topics", "Sldax", function(x) x@topics)

setGeneric("topics<-", function(x, value) standardGeneric("topics<-"))

#' @rdname Sldax
setMethod("topics<-", "Sldax", function(x, value) {
  x@topics <- value
  x
})

# theta
setGeneric("theta", function(x) standardGeneric("theta"))

#' @rdname Sldax
setMethod("theta", "Sldax", function(x) x@theta)

setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))

#' @rdname Sldax
setMethod("theta<-", "Sldax", function(x, value) {
  x@theta <- value
  x
})

# beta
setGeneric("beta_", function(x) standardGeneric("beta_"))

#' @rdname Sldax
setMethod("beta_", "Sldax", function(x) x@beta)

setGeneric("beta_<-", function(x, value) standardGeneric("beta_<-"))

#' @rdname Sldax
setMethod("beta_<-", "Sldax", function(x, value) {
  x@beta <- value
  x
})

# gamma
setGeneric("gamma_", function(x) standardGeneric("gamma_"))

#' @rdname Sldax
setMethod("gamma_", "Sldax", function(x) x@gamma)

setGeneric("gamma_<-", function(x, value) standardGeneric("gamma_<-"))

#' @rdname Sldax
setMethod("gamma_<-", "Sldax", function(x, value) {
  x@gamma <- value
  x
})

# alpha
setGeneric("alpha", function(x) standardGeneric("alpha"))

#' @rdname Sldax
setMethod("alpha", "Sldax", function(x) x@alpha)

setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))

#' @rdname Sldax
setMethod("alpha<-", "Sldax", function(x, value) {
  x@alpha <- value
  x
})

# ntopics
setGeneric("ntopics", function(x) standardGeneric("ntopics"))

#' @rdname Sldax
setMethod("ntopics", "Sldax", function(x) x@ntopics)

setGeneric("ntopics<-", function(x, value) standardGeneric("ntopics<-"))

#' @rdname Sldax
setMethod("ntopics<-", "Sldax", function(x, value) {
  x@ntopics <- value
  x
})

# nvocab
setGeneric("nvocab", function(x) standardGeneric("nvocab"))

#' @rdname Sldax
setMethod("nvocab", "Sldax", function(x) x@nvocab)

setGeneric("nvocab<-", function(x, value) standardGeneric("nvocab<-"))

#' @rdname Sldax
setMethod("nvocab<-", "Sldax", function(x, value) {
  x@nvocab <- value
  x
})

# Helper function (constructor) for Sldax class.
Sldax <- function(nvocab, topics, theta, beta, ntopics = 2.0, alpha = 1.0,
                  gamma = 1.0, ...) {
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

#' Summary functions for objects of class [Sldax-class]
#'
#' Obtain parameter estimates, model goodness-of-fit metrics, and posterior summaries.
#'
#' `get_zbar()` computes empirical topic proportions from slot `@topics`.
#' `est_theta()` estimates the mean or median theta matrix.
#' `est_beta()` estimates the mean or median beta matrix.
#' `get_toptopics()` creates a [`tibble`][tibble::tibble()] of the topic
#'    proportion estimates for the top `ntopics` topics per document sorted by
#'    probability.
#' `get_topwords()` creates a [`tibble`][tibble::tibble()] of topics and the
#'    top `nwords` words per topic sorted by probability or term score.
#' `get_coherence()` computes the coherence metric for each topic (see Mimno,
#'    Wallach, Talley, Leenders, & McCallum, 2011).
#' `get_exclusivity()` computes the exclusivity metric for each topic (see
#'    Roberts, Stewart, & Airoldi, 2013).
#' `post_regression()` creates a [`coda::mcmc`][coda::mcmc] object
#' containing posterior information for the regression model parameters.
#' `gg_coef()` plots regression coefficients (Warning: this function is
#'    deprecated. See `help("Deprecated")`.
#'
#' @name sldax-summary
NULL

#' @param mcmc_fit An object of class \linkS4class{Sldax}.
#' @param burn The number of draws to discard as a burn-in period (default: `0`).
#' @param thin The number of draws to skip as a thinning period (default: `1`; i.e., no thinning).
#' @param stat The summary statistic to use on the posterior draws (default: `"mean"`).
#'
#' @rdname sldax-summary
#' @export
setGeneric("est_beta",
           function(mcmc_fit,  burn = 0, thin = 1, stat = "mean") {

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if (length(dim(beta_(mcmc_fit))) != 3)
               stop("Only one draw of 'beta' available, so this function is not useful.")

             if ( !is.non_negative_integer(burn) ) stop("'burn' must be a non-negative integer.")
             if ( !is.positive_integer(thin) ) stop("'thin' must be a positive integer.")

             m <- nchain(mcmc_fit)
             if (burn >= m)
               stop("'burn' cannot exceed length of chain.")
             if (thin > (m - burn))
               stop("'thin' cannot exceed length of chain less 'burn'.")

             if (length(stat) > 1) {
               stat <- stat[1]
               message("Multiple arguments were supplied to 'stat'. Only using the first argument.")
             }
             if (!(stat %in% c("mean", "median")))
               stop("'stat' must be either 'mean' or 'median'.")

             standardGeneric("est_beta")
           }
)

#' @rdname sldax-summary
#' @export
setGeneric("est_theta",
           function(mcmc_fit, burn = 0, thin = 1, stat = "mean") {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if (length(dim(theta(mcmc_fit))) != 3)
               stop("Only one draw of 'theta' available, so this function is not useful.")

             if ( !is.non_negative_integer(burn) ) stop("'burn' must be a non-negative integer.")
             if ( !is.positive_integer(thin) ) stop("'thin' must be a positive integer.")

             m <- nchain(mcmc_fit)
             if (burn >= m)
               stop("'burn' cannot exceed length of chain.")
             if (thin > (m - burn))
               stop("'thin' cannot exceed length of chain less 'burn'.")

             if (length(stat) > 1) {
               stat <- stat[1]
               message("Multiple arguments were supplied to 'stat'. Only using the first argument.")
             }
             if (!(stat %in% c("mean", "median")))
               stop("'stat' must be either 'mean' or 'median'.")
             standardGeneric("est_theta")
           }
)

#' @param beta_ A \eqn{K} x \eqn{V} matrix of word-topic probabilities. Can be
#'   computed using [`est_beta()`][est_beta()]. Each row sums to 1.
#' @param docs The \eqn{D} x max(\eqn{N_d}) matrix of documents (word indices)
#'   used to fit the \linkS4class{Sldax} model.
#' @param nwords The number of highest-probability words per topic to consider where
#'   \eqn{M \le V} and \eqn{V} is the size of the corpus vocabulary. (default: `10`)
#'
#' @rdname sldax-summary
#' @export
setGeneric("get_coherence",
           function(beta_, docs, nwords = 10) {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(beta_))
               stop("Please supply an array to 'beta_'.")


             if ( length(dim(beta_)) != 2 )
               stop("'beta_' does not appear to be a K x V matrix.")

             if (any(beta_ < 0.0 | beta_ > 1.0)) stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             tol <- 0.001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if ( !is.positive_integer(nwords) ) stop("'nwords' must be a positive integer.")

             standardGeneric("get_coherence")
           }
)

#' @param weight The weight (between 0 and 1) to give to exclusivity (near 1) vs. frequency (near 0). (default: `0.7`)
#'
#' @rdname sldax-summary
#' @export
setGeneric("get_exclusivity",
           function(beta_, nwords = 10, weight = 0.7) {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(beta_))
               stop("Please supply an array to 'beta_'.")

             if (length(dim(beta_)) != 2)
               stop("'beta_' does not appear to be a K x V matrix.")

             if (any(beta_ < 0.0 | beta_ > 1.0)) stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             tol <- 0.001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if ( !is.positive_integer(nwords) ) stop("'nwords' must be a positive integer.")

             if ( ((weight >= 1.0) | (weight <= 0.0))) stop("'weight' must be between 0.0 and 1.0.")

             standardGeneric("get_exclusivity")
           }
)

#' @param theta A D x K matrix of K topic proportions for all D documents.
#' @param ntopics The number of topics to retrieve (default: all topics).
#'
#' @rdname sldax-summary
#' @export
setGeneric("get_toptopics",
           function(theta, ntopics) {

             if (missing(theta)) stop("Please supply a matrix to 'theta'.")
             if (!is.matrix(theta)) stop("Please supply a matrix to 'theta'.")
             if (any(theta < 0.0 | theta > 1.0)) stop("Entries of 'theta' must be between 0.0 and 1.0.")
             sum_rowsum_theta <- sum(rowSums(theta))
             d <- nrow(theta)
             K <- ncol(theta)
             tol <- 0.001
             if (sum_rowsum_theta > d + tol | sum_rowsum_theta < d - tol)
               stop("Rows of 'theta' must each sum to 1.0.")

             if (missing(ntopics)) ntopics <- K # Default
             if ( !is.positive_integer(ntopics) ) stop("'ntopics' must be a positive integer.")

             standardGeneric("get_toptopics")
           }
)

#' @param nwords The number of words to retrieve (default: all).
#' @param vocab A character vector of length V containing the vocabulary.
#' @param method If `"termscore"`, use term scores (similar to tf-idf). If
#'  `"prob"`, use probabilities (default: `"termscore"`).
#'
#' @return A \eqn{K} x \eqn{V} matrix of term-scores (comparable to tf-idf).
#' @rdname sldax-summary
#' @export
setGeneric("get_topwords",
           function(beta_, nwords, vocab, method = "termscore") {

             if (missing(beta_)) stop("Argument 'beta_' is missing.")
             if (!is.matrix(beta_)) stop("Please supply a matrix to 'beta_'.")
             if (any(beta_ < 0.0 | beta_ > 1.0))
               stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             V <- ncol(beta_)
             tol <- .001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if (missing(nwords)) nwords <- V # Default
             if ( !is.positive_integer(nwords) )
               stop("'nwords' must be a positive integer.")

             if (missing(vocab)) {
               vocab <- as.character(seq_len(V))
               message("'vocab' not supplied. Defaulting to indices 1, 2, ..., V.")
             }
             if (length(vocab) < 2L)
               stop("'vocab' must contain at least two elements.")
             if (length(vocab) != V)
               stop("The number of elements in 'vocab' should equal the number of columns in 'beta_'.")
             if (!is.character(vocab))
               stop("'vocab' must be a character vector.")
             if (nwords > length(unique(vocab)))
               stop("'nwords' cannot exceed the number of unique terms in 'vocab'.")

             if (length(method) > 1) {
               method <- method[1]
               message("Multiple arguments were supplied to 'method'. Only using the first argument.")
             }
             if (!(method %in% c("termscore", "prob")))
               stop("'method' must be either 'termscore' or 'prob'.")

             standardGeneric("get_topwords")
           }
)

#' @rdname sldax-summary
#' @export
setGeneric("get_zbar",
           function(mcmc_fit, burn = 0L, thin = 1L) {

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if ( !is.non_negative_integer(burn) )
               stop("'burn' must be a non-negative integer.")
             if ( !is.positive_integer(thin) ) stop("'thin' must be a positive integer.")

             m <- nchain(mcmc_fit)
             if (burn >= m)
               stop("'burn' cannot exceed length of chain.")
             if (thin > (m - burn))
               stop("'thin' cannot exceed length of chain less 'burn'.")

             standardGeneric("get_zbar")
           }
)

#' Summarize regression relationships for objects of class [Mlr-class]
#'
#' For SLDA or SLDAX models, label switching is handled during estimation in the
#' [`gibbs_sldax()`][gibbs_sldax()] function with argument `correct_ls`, so it
#' is not addressed by this function.
#'
#' @return An object of class [`coda::mcmc`][coda::mcmc] summarizing the posterior
#'   distribution of the regression coefficients and residual variance (if
#'   applicable). Convenience functions such as [`summary()`][coda::summary.mcmc()] and
#'   [`plot()`][coda::plot.mcmc()] can be used for posterior summarization.
#' @rdname sldax-summary
#' @export
setGeneric("post_regression",
           function(mcmc_fit) {

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if ( !( is(mcmc_fit, "Sldax") | is(mcmc_fit, "Mlr") | is(mcmc_fit, "Logistic")) )
               stop("'mcmc_fit' must be an `Sldax` or `Mlr` or `Logistic` object.")

             if ( extra(mcmc_fit)$call$model == "lda" )
               stop("The `lda` model does not contain regression parameters, so this function is not useful.")

             if ( nrow(eta(mcmc_fit)) < 2L )
               stop("Only one draw of 'eta' available, so this function is not useful.")

             standardGeneric("post_regression")
           }
)
