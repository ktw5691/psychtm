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
#' @param thin The period of iterations to keep after the burn-in period
#'   (default: 1).
#' @param docs A D x max(\eqn{N_d}) matrix of word indices for all documents.
#' @param V The number of unique terms in the vocabulary.
#' @param K The number of topics.
#' @param model A string denoting the type of model to fit. See 'Details'.
#'   (default: \code{"lda"}).
#' @param sample_beta A logical (default = \code{TRUE}): If \code{TRUE}, the
#'   topic-vocabulary distributions are sampled from their full conditional
#'   distribution.
#' @param sample_theta A logical (default = \code{TRUE}): If \code{TRUE}, the
#'   topic proportions will be sampled. CAUTION: This can be memory-intensive.
#' @param interaction_xcol EXPERIMENTAL: The column number of the design matrix for
#'   the additional predictors for which an interaction with the \eqn{K} topics is
#'   desired (default: \eqn{-1L}, no interaction). Currently only supports a single
#'   continuous predictor or a two-category categorical predictor represented as a
#'   single dummy-coded column.
#' @param alpha_ The hyper-parameter for the prior on the topic proportions
#'   (default: 1.0).
#' @param gamma_ The hyper-parameter for the prior on the topic-specific
#'   vocabulary probabilities (default: 1.0).
#' @param mu0 An optional q x 1 mean vector for the prior on the regression
#'   coefficients. See 'Details'.
#' @param sigma0 A q x q variance-covariance matrix for the prior on the
#'   regression coefficients. See 'Details'.
#' @param a0 The shape parameter for the prior on sigma2 (default: 0.001).
#' @param b0 The scale parameter for the prior on sigma2 (default: 0.001).
#' @param eta_start A q x 1 vector of starting values for the
#'   regression coefficients.
#' @param constrain_eta A logical (default = \code{FALSE}): If \code{TRUE}, the
#'   regression coefficients will be constrained so that they are in descending
#'   order; if \code{FALSE}, no constraints will be applied.
#' @param proposal_sd The proposal standard deviations for drawing the
#'   regression coefficients, N(0, proposal_sd(j)), \eqn{j = 1, \ldots, q}.
#'   Only used for \code{model = "slda_logit"} and
#'   \code{model = "sldax_logit"} (default: 2.38 for all coefficients).
#' @param return_assignments A logical (default = \code{FALSE}): If
#'   \code{TRUE}, returns an N x \eqn{max N_d} x M array of topic assignments
#'   in slot @topics. CAUTION: this can be memory-intensive.
#' @param verbose Should parameter draws be output during sampling? (default:
#'   \code{FALSE}).
#' @param display_progress Should percent progress of sampler be displayed
#'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
#'   \code{display_progress} be set to \code{TRUE} at any given time.
#'
#' @return An object of class \code{\linkS4class{Sldax}}.
#' @family Gibbs sampler
gibbs_sldax <- function(formula, data, m = 100, burn = 0, thin = 1,
                        docs, V, K = 2L,
                        model = c("lda", "slda", "sldax",
                                  "slda_logit", "sldax_logit"),
                        sample_beta = TRUE, sample_theta = TRUE,
                        interaction_xcol = -1L, alpha_ = 1.0, gamma_ = 1.0,
                        mu0 = NULL, sigma0 = NULL, a0 = NULL, b0 = NULL,
                        eta_start = NULL, constrain_eta = FALSE,
                        proposal_sd = NULL, return_assignments = FALSE,
                        correct_ls = TRUE,
                        verbose = FALSE, display_progress = FALSE) {

  # Start timing
  t1 <- Sys.time()

  # Check if any arguments without defaults were not supplied by user
  if (missing(formula) & model != "lda") missing_msg(formula)
  if (missing(data) & model != "lda") missing_msg(data)
  if (missing(docs)) missing_msg(docs)
  if (missing(V)) missing_msg(V)

  # Check for non-negative integers only in docs
  if (any(!is.non_negative_integer(docs)))
    stop("Invalid elements in 'docs'. 'docs' can only contain non-negative integers")
  # Check for documents of at least two words in length in docs
  if (any(apply(docs, 1, function(x) length(x[x > 0])) < 2))
    stop("Each document (row) in 'docs' must contain at least 2 positive integers")
  # Check for missing values in documents
  if (any(is.na(docs))) stop("NA found in 'docs'. Please use 0 instead to indicate unused word positions")

  # Check that data and docs have same number of observations
  if (model != "lda") {
    if (NROW(data) != NROW(docs))
      stop("'data' and 'docs' have unequal numbers of observations")
    # Check for missing values in data
    if (any(is.na(data))) stop("Cannot handle missing values in 'data'")
  }

  # Check model and convert to integer codes
  if (is.character(model)) {
    if (length(model) == 1) {
      model <- model
    } else {
      model <- model[1]
    }
    if (!model %in% c("lda", "slda", "sldax", "slda_logit", "sldax_logit")) {
      print(model)
      stop("'model' not recognized")
    }
  } else {
    print(model)
    stop("'model' not recognized")
  }
  model <- switch(model,
                  lda = 1L,
                  slda = 2L,
                  sldax = 3L,
                  slda_logit = 4L,
                  sldax_logit = 5L)

  sldax_call <- match.call()

  if (model != 1L) {
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
        x <- x[, -ip]
      }
      x <- as.matrix(x) # If y ~ 1 supplied, x has dims D x 0
      if (dim(x)[2] == 0 & model %in% c(3, 5))
        stop("Design matrix 'x' has 0 columns.
              Possible reason: don't supply y ~ 1 or y ~ 0 or y ~ -1 as
              'formula'.")

      if (model == 4 | model == 5) { # y should be dichotomous
        # Check if y is factor and convert to 0/1 numeric if it is
        if (is.factor(y)) y <- as.numeric(y) - 1
        # Check that y is 0/1
        if (length(unique(y)) == 2) {
          if (sum(unique(y) %in% c(0, 1)) != 2)
            stop("'y' has 2 unique values, but they cannot be coerced to 0/1")
        } else {
          stop("'y' has > 2 unique values and is not dichotomous")
        }
        y <- as.matrix(y)
        attr(y, "names") <- NULL
        if (sum(unique(y) %in% c(0, 1)) != 2)
          stop("'y' has 2 unique values, but they cannot be coerced to 0/1")
      } else {
        y <- as.matrix(y) # Assume y is continuous
      }
    }
  }

  if (model == 1) {
    if (!is.null(mu0) | !is.null(sigma0) | !is.null(a0) | !is.null(b0) |
        !is.null(eta_start) | !is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'lda' model")
    }
  }
  if (model %in% c(2, 3)) {
    if (!is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'slda'/'sldax' model")
    }
  }
  if (model %in% c(4, 5)) {
    if (!is.null(a0) | !is.null(b0)) {
      stop("Invalid parameter supplied for 'slda_logit'/'sldax_logit' model")
    }
  }

  # Default priors for supervised models
  if (model %in% c(2, 4)) { # slda or slda_logit
    p <- 0
    q_ <- K
    if (is.null(mu0)) mu0 <- rep(0, q_)
    if (is.null(sigma0)) {
      if (model == 2) sigma0 <- diag(1, q_)
      if (model == 4) sigma0 <- diag(2.5 ^ 2, q_)
    }
    if (is.null(eta_start)) eta_start <- seq(2, -2, length.out = q_)
  }
  if (model %in% c(3, 5)) { # sldax or sldax_logit
    if (interaction_xcol <= 0) {
      p <- NCOL(x)
    } else {
      p <- ncol(x) + K - 1
    }
    q_ <- p + K
    if (is.null(mu0)) mu0 <- rep(0, q_)
    if (is.null(sigma0)) {
      if (model == 3) sigma0 <- diag(1, q_)
      if (model == 5) sigma0 <- diag(2.5 ^ 2, q_)
    }
    if (is.null(eta_start))
      eta_start <- c(rep(0, ncol(x)), seq(2, -2, length.out = q_ - ncol(x)))
  }
  if (model %in% c(2, 3)) {
    if (is.null(a0)) a0 <- 0.001
    if (is.null(b0)) b0 <- 0.001
  }
  if (model == 4) {
    q_ <- K
    if (is.null(proposal_sd)) proposal_sd <- rep(2.38, q_)
  }
  if (model == 5) {
    if (interaction_xcol <= 0) {
      p <- ncol(x)
    } else {
      p <- ncol(x) + K - 1
    }
    q_ <- p + K
    if (is.null(proposal_sd)) proposal_sd <- rep(2.38, q_)
  }

  # Set up unused defaults for C++ call
  if (model == 1) {
    y           <- -1L
    x           <- matrix(0, 1, 1)
    mu0         <- 0
    sigma0      <- matrix(0, 1, 1)
    a0          <- -1
    b0          <- -1
    eta_start   <- 0
    proposal_sd <- 0
  }
  if (model == 2) {
    x           <- matrix(0, 1, 1)
    proposal_sd <- 0
  }
  if (model == 3) {
    proposal_sd <- 0
  }
  if (model == 4) {
    x  <- matrix(0, 1, 1)
    a0 <- -1
    b0 <- -1
  }
  if (model == 5) {
    a0 <- -1
    b0 <- -1
  }

  # Check m
  chk_m <- is.positive_integer(m)
  if (!chk_m) stop("'m' is not a positive integer")

  # Check burn
  chk_burn <- is.non_negative_integer(burn)
  if (!chk_burn) stop("'burn' is not a non-negative integer")

  # Check number of topics
  chk_K <- is.positive_integer(K)
  if (!chk_K) stop("'K' is not a positive integer")

  # Check sample_beta
  chk_sample_beta <- check_logical(sample_beta)
  if (!chk_sample_beta) stop("'sample_beta' is not TRUE/FALSE")

  # Check sample_theta
  chk_sample_theta <- check_logical(sample_theta)
  if (!chk_sample_theta) stop("'sample_theta' is not TRUE/FALSE")

  # Check constrain_eta
  chk_constrain <- check_logical(constrain_eta)
  if (!chk_constrain) stop("'constrain_eta' is not TRUE/FALSE")

  # Check verbose
  chk_verbose <- check_logical(verbose)
  if (!chk_verbose) stop("'verbose' is not an TRUE/FALSE")

  # Check display_progress
  chk_display <- check_logical(display_progress)
  if (!chk_display) stop("'display_progress' is not TRUE/FALSE")

  # Check correct_ls
  chk_correct_ls <- check_logical(correct_ls)
  if (!chk_correct_ls) stop("'correct_ls' is not TRUE/FALSE")

  # Check for requirements to correct label switching
  if (isTRUE(correct_ls) & !(isTRUE(sample_beta) & isTRUE(sample_theta)))
    stop("'sample_beta' and 'sample_theta' must both be TRUE to correct label switching")

  # TODO: Currently does not handle label switching correction when topic-covariate interaction specified
  if (isTRUE(correct_ls) & isFALSE(interaction_xcol < 1))
    stop("Label switching correction is not currently supported for models with topic-covariate interactions")

  res_out <- tryCatch({
    res <- .gibbs_sldax_cpp(docs = docs, V = V, m = m, burn = burn, thin = thin,
                            K = K, model = model, y = y, x = x,
                            mu0 = mu0, sigma0 = sigma0, a0 = a0, b0 = b0,
                            eta_start = eta_start, proposal_sd = proposal_sd,
                            interaction_xcol = interaction_xcol,
                            alpha_ = alpha_, gamma_ = gamma_,
                            constrain_eta = constrain_eta,
                            sample_beta = sample_beta,
                            sample_theta = sample_theta,
                            return_assignments = return_assignments,
                            verbose = verbose,
                            display_progress = display_progress)
    if (model %in% c(3, 5) & !is.null(colnames(x))) {
      if (interaction_xcol < 1) {
        colnames(eta(res)) <- c(colnames(x), paste0("topic", seq_len(K)))
      }
      else {
        colnames(eta(res)) <- c(colnames(x), paste0("topic", seq_len(K)),
                                paste0(colnames(x)[interaction_xcol], ":topic",  seq_len(K - 1)))
      }
    } else if (model %in% c(2, 4)) {
      colnames(eta(res)) <- paste0("topic", seq_len(K))
    }

    # Correct label switching; requires sampling theta and beta
    #   TODO: Label switching not implemented for models with topic-covariate interactions
    correct_label_switch <- FALSE
    if (correct_ls) {
      # relabel_out$permutations is Nchain x K array of permutation indices
      relabel_out = stephens(aperm(theta(res), c(3, 1, 2)), maxiter = 100) # First arg needs to be Nchain x D x K array
      if (relabel_out$iterations >= 100) {
        cat("Relabeling failed to converge, no label switching correction applied.")
      } else {
        # Permute theta
        perm_theta <- permute.mcmc(aperm(theta(res), c(3, 2, 1)), relabel_out$permutations) # First arg needs to be Nchain x K x npar array; npar won't be swapped, but rows (K) will be permuted
        perm_theta <- aperm(perm_theta$output, c(3, 2, 1)) # Return to D x K x Nchain format
        theta(res) <- perm_theta
        rm(perm_theta)

        # Permute beta
        perm_beta <- permute.mcmc(aperm(beta_(res), c(3, 1, 2)), relabel_out$permutations)
        perm_beta <- aperm(perm_beta$output, c(2, 3, 1)) # Return to K x V x Nchain format
        beta_(res) <- perm_beta
        rm(perm_beta)

        if (model != 1) { # Skip for LDA model
          # Permute eta for topics
          perm_eta <- permute.mcmc(array(eta(res)[, seq(p + 1, q_)], dim = c(nchain(res), K, 1)), relabel_out$permutations)
          perm_eta <- array(perm_eta$output, c(nchain(res), K))
          eta(res)[, seq(p + 1, q_)] <- perm_eta
          rm(perm_eta)
        }

        # Permute topic assignments for all words
        if (return_assignments) {
          perm_z <- topics(res)
          for (iter in seq_len(nchain(res))){
            perm_i <- relabel_out$permutations[iter, ]
            if (isTRUE(all.equal(seq_len(K), perm_i))) next # No change needed, go to next sample
            old_z <- perm_z[, , iter] # Original draws to check
            for (k in seq_len(K)) {
              new_index = perm_i[k]
              if (new_index == k) next # Index does not need to change for this topic, go to next topic
              perm_z[, , iter][old_z == k] <- new_index
            }
          }
          topics(res) <- perm_z
        }
        correct_label_switch <- TRUE
      }
    }
    t2 <- Sys.time()
    extra(res) <- list(time_elapsed = t2 - t1,
                       start_time = t1,
                       end_time = t2,
                       corrected_label_switching = correct_label_switch,
                       call = sldax_call)
    return(res)
  }, error = function(err_cond) {
    stop(err_cond)
  })

  return(res_out)
}

#' Fit linear regression model
#'
#' \code{gibbs_mlr} is used to fit a Bayesian linear regression model using
#' Gibbs sampling.
#'
#' For \code{mu0}, by default, we use a vector of \eqn{p} 0s for \eqn{p}
#' regression coefficients.
#'
#' For \code{sigma0}, by default, we use a \eqn{p} x \eqn{p} identity matrix.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted.
#' @param data An optional data frame containing the variables in the model.
#' @param m The number of iterations to run the Gibbs sampler (default: 100).
#' @param burn The number of iterations to discard as the burn-in period
#'   (default: 0).
#' @param thin The period of iterations to keep after the burn-in period
#'   (default: 1).
#' @param mu0 An optional p x 1 mean vector for the prior on the regression
#'   coefficients. See 'Details'.
#' @param sigma0 A p x p variance-covariance matrix for the prior on the
#'   regression coefficients. See 'Details'.
#' @param a0 The shape parameter for the prior on sigma2 (default: 0.001).
#' @param b0 The scale parameter for the prior on sigma2 (default: 0.001).
#' @param eta_start A p x 1 vector of starting values for the regression
#'   coefficients.
#' @param verbose Should parameter draws be output during sampling? (default:
#'   \code{FALSE}).
#' @param display_progress Should percent progress of sampler be displayed
#'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
#'   \code{display_progress} be set to \code{TRUE} at any given time.
#'
#' @return An object of class \code{\linkS4class{Mlr}}.
#' @family Gibbs sampler
gibbs_mlr <- function(formula, data, m = 100, burn = 0, thin = 1,
                      mu0 = NULL, sigma0 = NULL, a0 = NULL, b0 = NULL,
                      eta_start = NULL,
                      verbose = FALSE, display_progress = FALSE) {

  # Start timing
  t1 <- Sys.time()

  # Check if any arguments without defaults were not supplied by user
  if (missing(formula)) missing_msg(formula)
  if (missing(data)) missing_msg(data)

  # Check for missing values in data
  if (any(is.na(data))) stop("Cannot handle missing values in 'data'")

  mlr_call <- match.call()
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
    y <- as.matrix(y)
  }

  # Default priors
  q_ <- ncol(x)
  if (is.null(mu0)) mu0 <- rep(0, q_)
  if (is.null(sigma0)) sigma0 <- diag(1, q_)
  if (is.null(eta_start)) eta_start <- rep(0, q_)
  if (is.null(a0)) a0 <- 0.001
  if (is.null(b0)) b0 <- 0.001

  # Check m
  chk_m <- is.positive_integer(m)
  if (!chk_m) stop("'m' is not a positive integer")

  # Check burn
  chk_burn <- is.non_negative_integer(burn)
  if (!chk_burn) stop("'burn' is not a non-negative integer")

  # Check verbose
  chk_verbose <- check_logical(verbose)
  if (!chk_verbose) stop("'verbose' is not an TRUE/FALSE")

  # Check display_progress
  chk_display <- check_logical(display_progress)
  if (!chk_display) stop("'display_progress' is not TRUE/FALSE")

  res_out <- tryCatch({
    res <- .gibbs_mlr_cpp(m, burn, thin, y, x, mu0, sigma0, eta_start, a0, b0,
                          verbose, display_progress)
    colnames(eta(res)) <- colnames(x)
    t2 <- Sys.time()
    extra(res) <- list(time_elapsed = t2 - t1,
                       start_time = t1,
                       end_time = t2,
                       call = mlr_call)
    return(res)
    }, error = function(err_cond) {
      stop(err_cond)
    })

  return(res_out)
}

#' Fit logistic regression model
#'
#' \code{gibbs_logistic} is used to fit a Bayesian logistic regression model
#' using Gibbs sampling.
#'
#' For \code{mu0}, by default, we use a vector of \eqn{p} 0s for \eqn{p}
#' regression coefficients.
#'
#' For \code{sigma0}, by default, we use a \eqn{p} x \eqn{p} diagonal matrix
#' with diagonal elements (variances) of 6.25.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted.
#' @param data An optional data frame containing the variables in the model.
#' @param m The number of iterations to run the Gibbs sampler (default: 100).
#' @param burn The number of iterations to discard as the burn-in period
#'   (default: 0).
#' @param thin The period of iterations to keep after the burn-in period
#'   (default: 1).
#' @param mu0 An optional p x 1 mean vector for the prior on the regression
#'   coefficients. See 'Details'.
#' @param sigma0 A p x p variance-covariance matrix for the prior on the
#'   regression coefficients. See 'Details'.
#' @param eta_start A p x 1 vector of starting values for the
#'   regression coefficients.
#' @param proposal_sd The proposal standard deviations for drawing the
#'   regression coefficients, N(0, proposal_sd(j)), \eqn{j = 1, \ldots, p}
#'   (default: 2.38 for all coefficients).
#' @param verbose Should parameter draws be output during sampling? (default:
#'   \code{FALSE}).
#' @param display_progress Should percent progress of sampler be displayed
#'   (default: \code{FALSE}). Recommended that only one of \code{verbose} and
#'   \code{display_progress} be set to \code{TRUE} at any given time.
#'
#' @return An object of class \code{\linkS4class{Logistic}}.
#' @family Gibbs sampler
gibbs_logistic <- function(formula, data, m = 100, burn = 0, thin = 1,
                           mu0 = NULL, sigma0 = NULL,
                           eta_start = NULL, proposal_sd = NULL,
                           verbose = FALSE, display_progress = FALSE) {

  # Start timing
  t1 <- Sys.time()

  # Check if any arguments without defaults were not supplied by user
  if (missing(formula)) missing_msg(formula)
  if (missing(data)) missing_msg(data)

  # Check for missing values in data
  if (any(is.na(data))) stop("Cannot handle missing values in 'data'")

  logistic_call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mind <- match(c("formula", "data"), names(mf), 0L)
  if (sum(mind) > 0) {
    mf <- mf[c(1L, mind)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf, "any")
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)

    # Check if y is factor and convert to 0/1 numeric if it is
    if (is.factor(y)) y <- as.numeric(y) - 1
    # Check that y is 0/1
    if (length(unique(y)) == 2) {
      if (sum(unique(y) %in% c(0, 1)) != 2)
        stop("'y' has 2 unique values, but they cannot be coerced to 0/1")
    } else {
      stop("'y' has > 2 unique values and is not dichotomous")
    }
    y <- as.matrix(y)
    attr(y, "names") <- NULL
    if (sum(unique(y) %in% c(0, 1)) != 2)
      stop("'y' has 2 unique values, but they cannot be coerced to 0/1")
  }

  # Default priors
  q_ <- ncol(x)
  if (is.null(mu0)) mu0 <- rep(0, q_)
  if (is.null(sigma0)) sigma0 <- diag(2.5 ^ 2, q_)
  if (is.null(eta_start)) eta_start <- rep(0, q_)
  if (is.null(proposal_sd)) proposal_sd <- rep(2.38, q_)

  # Check m
  chk_m <- is.positive_integer(m)
  if (!chk_m) stop("'m' is not a positive integer")

  # Check burn
  chk_burn <- is.non_negative_integer(burn)
  if (!chk_burn) stop("'burn' is not a non-negative integer")

  # Check verbose
  chk_verbose <- check_logical(verbose)
  if (!chk_verbose) stop("'verbose' is not an TRUE/FALSE")

  # Check display_progress
  chk_display <- check_logical(display_progress)
  if (!chk_display) stop("'display_progress' is not TRUE/FALSE")

  # Check proposal_sd
  chk_proposal_sd <- (length(proposal_sd) == q_) &
    (sum(sign(proposal_sd) > 0) == q_)
  if (!chk_proposal_sd)
    stop(paste("'proposal_sd' must be a positive vector of length", q_))

  res_out <- tryCatch({
    res <- .gibbs_logistic_cpp(m, burn, thin, y, x, mu0, sigma0,
                               eta_start, proposal_sd,
                               verbose, display_progress)
    colnames(eta(res)) <- colnames(x)
    t2 <- Sys.time()
    extra(res) <- list(time_elapsed = t2 - t1,
                       start_time = t1,
                       end_time = t2,
                       call = logistic_call)
    return(res)
  }, error = function(err_cond) {
    stop(err_cond)
  })

  return(res_out)
}
