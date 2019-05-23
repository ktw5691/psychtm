#' Fit supervised or unsupervised topic models (sLDAX or LDA)
#'
#' \code{gibbs_sldax} is used to fit both supervised and unsupervised topic models.
#'
#' The number of regression coefficients q in supervised topic models is
#' determined as follows: For the sLDA model with only the \eqn{K} topics as
#' predictors, \eqn{q = K}; for the sLDAX model with \eqn{K} topics and \eqn{p}
#' additional predictors, there are two possibilities: (1) If no interaction
#' between an additional covariate and the \eqn{K} topics is desired
#' (default: \code{interaction_xcol = -1L}), \eqn{q = p + K}; (2) if an
#' interaction between an additional covariate and the \eqn{K} topics is desired
#' (e.g., \code{interaction_xcol = 1}), \eqn{q = p + 2K - 1}. If you supply
#' custom values for prior parameters \code{mu0} or \code{sigma0}, be sure that
#' the length of \code{mu0} (\eqn{q}) and/or the number of rows and columns of
#' \code{sigma0} (\eqn{q \times q}) are correct. If you supply custom starting
#' values for \code{eta_start}, be sure that the length of \code{eta_start} is
#' correct.
#'
#' For \code{model}, one of \code{c("lda", "slda", "sldax", "slda_logit",
#' "sldax_logit")} are possible. \code{"lda"}: unsupervised topic model;
#' \code{"slda"}: supervised topic model with a continuous outcome;
#' \code{"sldax"}: supervised topic model with a continuous outcome and
#' additional predictors of the outcome; \code{"slda_logit"}: supervised topic
#' model with a dichotomous outcome (0/1); \code{"sldax_logit"}: supervised
#' topic model with a dichotomous outcome (0/1) and additional predictors of the
#' outcome.
#'
#' For \code{mu0}, the first \eqn{p} elements correspond to coefficients for the
#' \eqn{p} additional predictors (if none, \eqn{p = 0}), while elements
#' \eqn{p + 1} to \eqn{p + K} correspond to coefficients for the \eqn{K} topics,
#' and elements \eqn{p + K + 1} to \eqn{p + 2K - 1} correspond to coefficients
#' for the interaction (if any) between one additional predictor and the \eqn{K}
#' topics. By default, we use a vector of \eqn{q} 0s.
#'
#' For \code{sigma0}, the first \eqn{p} rows/columns correspond to coefficients
#' for the \eqn{p} additional predictors (if none, \eqn{p = 0}), while
#' rows/columns \eqn{p + 1} to \eqn{p + K} correspond to coefficients for the
#' \eqn{K} topics, and rows/columns \eqn{p + K + 1} to \eqn{p + 2K - 1}
#' correspond to coefficients for the interaction (if any) between one
#' additional predictor and the \eqn{K} topics. By default, we use an identity
#' matrix for \code{model = "slda"} and \code{model = "sldax"} and a diagonal
#' matrix with diagonal elements (variances) of 6.25 for
#' \code{model = "slda_logit"} and \code{model = "sldax_logit"}.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted.
#' @param data An optional data frame containing the variables in the model.
#' @param m The number of iterations to run the Gibbs sampler (default: 100).
#' @param burn The number of iterations to discard as the burn-in period
#'   (default: 0).
#' @param docs A D x max(\eqn{N_d}) matrix of word indices for all documents.
#' @param w A D x V matrix of counts for all documents and vocabulary terms.
#' @param K The number of topics.
#' @param model A string denoting the type of model to fit. See 'Details'.
#'   (default: \code{"lda"}).
#' @param y An optional D x 1 vector of binary outcomes (0/1) to be predicted.
#'   Not needed if \code{data} is supplied.
#' @param x An optional D x p design matrix of additional predictors (do NOT
#'   include an intercept column!). Not needed if \code{data} is supplied.
#' @param interaction_xcol The column number of the design matrix for the
#' additional predictors for which an interaction with the \eqn{K} topics is desired
#' (default: \eqn{-1L}, no interaction). Currently only supports a single continuous
#' predictor or a two-category categorical predictor represented as a single
#' dummy-coded column.
#' @param alpha_ The hyper-parameter for the prior on the topic proportions
#'   (default: 0.1).
#' @param gamma_ The hyper-parameter for the prior on the topic-specific
#'   vocabulary probabilities (default: 1.01).
#' @param mu0 An optional q x 1 mean vector for the prior on the regression
#'   coefficients. See 'Details'.
#' @param sigma0 A q x q variance-covariance matrix for the prior on the
#'   regression coefficients. See 'Details'.
#' @param eta_start A q x 1 vector of starting values for the
#'   regression coefficients.
#' @param constrain_eta A logical (default = \code{TRUE}): If \code{TRUE}, the
#'   regression coefficients will be constrained so that they are in descending
#'   order; if \code{FALSE}, no constraints will be applied.
#' @param proposal_sd The proposal standard deviations for drawing the
#'   regression coefficients, N(0, proposal_sd(j)), \eqn{j = 1, \ldots, q}.
#'   Only used for \code{model = "slda_logit"} and
#'   \code{model = "sldax_logit"} (default: 2.38 for all coefficients).
#' @param verbose Should parameter draws be output during sampling? (default:
#'   \code{FALSE}).
#' @param display_progress Should percent progress of sampler be displayed
#'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
#'   \code{display_progress} be set to \code{TRUE} at any given time.
#'
#' @return An object of class \code{\linkS4class{Sldax}}.
gibbs_sldax = function(formula, data, m = 100, burn = 0, docs, w, K = 2L,
                       model = c("lda", "slda", "sldax", "slda_logit", "sldax_logit"),
                       y = NULL, x = NULL, interaction_xcol = -1L,
                       alpha_ = 0.1, gamma_ = 1.01,
                       mu0 = NULL, sigma0 = NULL, a0 = NULL, b0 = NULL,
                       eta_start = NULL, constrain_eta = TRUE,
                       proposal_sd = NULL,
                       verbose = FALSE, display_progress = FALSE) {

  mf <- match.call(expand.dots = FALSE)
  mind <- match(c("formula", "data"), names(mf), 0L)
  if (sum(mind) > 0) {
    mf <- mf[c(1L, mind)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf, "numeric")
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    if ("(Intercept)" %in% dimnames(x)[[2]]) {
      ip <- match("(Intercept)", dimnames(x)[[2]])
      x = x[, -1]
    }
    x = as.matrix(x) # If y ~ 1 supplied, x has dims D x 0
    if (dim(x)[2] == 0)
      stop("Design matrix `x` has 0 columns. Possible reason: don't supply y ~ 1 or y ~ 0 or y ~ -1 as `formula`.")
    y = as.matrix(y)
  }

  # Check model and convert to integer codes
  if (is.character(model)) {
    if (length(model) == 1) {
      model = model
    } else {
      model = model[1]
    }
    if (!model %in% c("lda", "slda", "sldax", "slda_logit", "sldax_logit")) {
      print(model)
      stop("'model' not recognized")
    }
  } else {
    print(model)
    stop("'model' not recognized")
  }
  model = switch(model,
                 lda = 1L,
                 slda = 2L,
                 sldax = 3L,
                 slda_logit = 4L,
                 sldax_logit = 5L)

  if (model == 1) {
    if (!is.null(y) || !is.null(x) || !is.null(mu0) || !is.null(sigma0) ||
        !is.null(a0) || !is.null(b0) || !is.null(eta_start) ||
        !is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'lda' model")
    }
  }
  if (model == 2) {
    if (!is.null(x) || !is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'slda' model")
    }
  }
  if (model == 3) {
    if (!is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'sldax' model")
    }
  }
  if (model == 4) {
    if (!is.null(x) || !is.null(a0) || !is.null(b0)) {
      stop("Invalid parameter supplied for 'slda_logit' model")
    }
  }
  if (model == 5) {
    if (!is.null(a0) || !is.null(b0)) {
      stop("Invalid parameter supplied for 'sldax_logit' model")
    }
  }

  # Default priors for supervised models
  if (model %in% c(2, 4)) {
    q_ = K
    if (is.null(mu0)) mu0 = rep(0, q_)
    if (is.null(sigma0)) {
      if (model == 2) sigma0 = diag(1, q_)
      if (model == 4) sigma0 = diag(2.5 ^ 2, q_) # Roughly equivalent to recommended Cauchy prior scales of 2.5 from Gelman paper
    }
    if (is.null(eta_start)) eta_start = seq(2, -2, length.out = q_)
  }
  if (model %in% c(3, 5)) {
    if (interaction_xcol <= 0) {
      p = ncol(x)
    } else {
      p = ncol(x) + K - 1
    }
    q_ = p + K
    if (is.null(mu0)) mu0 = rep(0, q_)
    if (is.null(sigma0)) {
      if (model == 3) sigma0 = diag(1, q_)
      if (model == 5) sigma0 = diag(2.5 ^ 2, q_) # Roughly equivalent to recommended Cauchy prior scales of 2.5 from Gelman paper
    }
    if (is.null(eta_start)) eta_start = c(rep(0, ncol(x)), seq(2, -2, length.out = q_ - ncol(x)))
  }
  if (model %in% c(2, 3)) {
    if (is.null(a0)) a0 = 0.001
    if (is.null(b0)) b0 = 0.001
  }
  if (model == 4) {
    q_ = K
    if (is.null(proposal_sd)) proposal_sd = rep(2.38, q_)
  }
  if (model == 5) {
    if (interaction_xcol <= 0) {
      p = ncol(x)
    } else {
      p = ncol(x) + K - 1
    }
    q_ = p + K
    if (is.null(proposal_sd)) proposal_sd = rep(2.38, q_)
  }

  # Set up unused defaults for C++ call
  if (model == 1) {
    y = -1L
    x = matrix(0, 1, 1)
    mu0 = 0
    sigma0 = matrix(0, 1, 1)
    a0 = -1
    b0 = -1
    eta_start = 0
    proposal_sd = 0
  }
  if (model == 2) {
    x = matrix(0, 1, 1)
    proposal_sd = 0
  }
  if (model == 3) {
    proposal_sd = 0
  }
  if (model == 4) {
    x = matrix(0, 1, 1)
    a0 = -1
    b0 = -1
  }
  if (model == 5) {
    a0 = -1
    b0 = -1
  }

  # Check m
  if (is.numeric(m)) {
    if (!is.na(m) & !is.infinite(m)) {
      if (m %% 1 == 0) {
        m = as.integer(m)
      }
    } else {
      print(m)
      stop("'m' is not an integer")
    }
  } else {
    print(m)
    stop("'m' is not an integer")
  }

  # Check burn
  if (is.numeric(burn)) {
    if (!is.na(burn) & !is.infinite(burn)) {
      if (burn %% 1 == 0) {
        burn = as.integer(burn)
      }
    } else {
      print(burn)
      stop("'burn' is not an integer")
    }
  } else {
    print(burn)
    stop("'burn' is not an integer")
  }

  # Check number of topics
  if (is.numeric(K)) {
    if (!is.na(K) & !is.infinite(K)) {
      if (K %% 1 == 0) {
        K = as.integer(K)
      }
    } else {
      print(K)
      stop("'K' is not an integer")
    }
  } else {
    print(K)
    stop("'K' is not an integer")
  }

  gibbs_sldax_cpp(docs, w, m, burn, K, model, y, x, mu0, sigma0, a0, b0,
                  eta_start, proposal_sd, interaction_xcol, alpha_, gamma_,
                  constrain_eta, verbose, display_progress)
}
