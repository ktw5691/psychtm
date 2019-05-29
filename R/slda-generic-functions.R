#' Generic function to retrieve the top topics for each document
#'
#' Functions to compute and visualize the most probable topics in a document or
#' the most probable words in a topic.
#'
#' @param mcmc_fit An Lda object.
#' @param ntopics The number of topics to retrieve.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param stat The summary statistic to use on the posterior draws. Default: \code{"mean"}.
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_toptopics",
           function(mcmc_fit, ntopics, burn = 0, thin = 1, stat = "mean") {

             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")
             if (!isClass(mcmc_fit, "Lda"))
               stop("'mcmc_fit' must be an Lda object.")
             if ( (burn %% 1) != 0) stop("'burn' must be an integer.")
             if (burn < 0) stop("'burn' must be non-negative.")
             if ( (thin %% 1) != 0) stop("'thin' must be an integer.")
             if (thin < 1) stop("'thin' must be positive.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("'burn' cannot exceed length of chain.")
             if (thin >= (m - burn))
               stop("'thin' cannot exceed length of chain less 'burn'.")

             if (length(stat > 1)) stat <- stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("stat must be either 'mean' or 'median'")
             if (is.null(stat)) stat <- "mean" # Default to mean

             standardGeneric("get_toptopics")
           }
)

#' Generic function to retrieve the most probable words for each topic
#'
#' @param beta_ A \eqn{K} x \eqn{V} matrix of word-topic probabilities where
#'   each row sums to 1.
#' @param nwords The number of words to retrieve.
#' @param vocab A character vector containing the vocabulary.
#' @param method If "termscore", use term scores (similar to tf-idf). If "prob",
#'   use probabilities. Default: "termscore".
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_topwords",
           function(beta_, nwords, vocab, method = "termscore") {

             if ( (nwords %% 1) != 0) stop("'nwords' must be an integer.")
             if (nwords < 1) stop("'nwords' must be an integer greater than 0.")
             if (length(vocab) == 0)
               stop("'vocab' must contain at least one element.")
             if (nwords > length(unique(vocab)))
               stop("'n_words' cannot exceed the number of unique terms in
                     'vocab'.")

             if (length(method > 1)) method <- method[1]
             if (!(method %in% c("termscore", "prob")))
               stop("'stat' must be either 'termscore' or 'prob'")
             if (is.null(method)) method <- "termscore" # Default to termscore

             standardGeneric("get_topwords")
           }
)

#' Generic function to retrieve the empirical topic proportions
#'
#' Compute empirical topic proportions (zbar) from \code{@topics} in a
#' \code{Lda} object.
#' @rdname slda-gettop-methods
setGeneric("get_zbar",
           function(mcmc_fit, burn = 0, thin = 1) {

             if (!isClass(mcmc_fit, "Lda"))
               stop("'mcmc_fit' must be an Lda object.")
             if (is.null(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")

             if ( (burn %% 1) != 0) stop("'burn' must be an integer.")
             if (burn < 0) stop("'burn' must be non-negative.")
             if ( (thin %% 1) != 0) stop("'thin' must be an integer.")
             if (thin < 1) stop("'thin' must be positive.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("'burn' cannot exceed length of chain.")
             if (thin >= (m - burn))
               stop("'thin' cannot exceed length of chain less 'burn'.")

             standardGeneric("get_zbar")
           }
)

#' Generic function to plot the regression coefficients for sLDAX models
#'
#' @param varnames A character vector of variable names for additional
#'   predictors (if any).
#' @param errorbw Controls the width of the +/- 2 posterior standard error bars.
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("gg_coef",
           function(mcmc_fit, beta_, nwords, vocab, varnames, burn = 0,
                    thin = 1, method = "termscore", stat = "mean",
                    errorbw = 0.5) {

             if (!isClass(mcmc_fit, "Lda"))
               stop("'mcmc_fit' must be an Lda object.")
             if (is.null(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if ( (nwords %% 1) != 0) stop("'nwords' must be an integer.")
             if (nwords < 1) stop("'nwords' must be an integer greater than 0.")
             if (length(vocab) == 0)
               stop("'vocab' must contain at least one element.")
             if (nwords > length(unique(vocab)))
               stop("'n_words' cannot exceed the number of unique terms in
                     'vocab'.")
             if ( (burn %% 1) != 0) stop("'burn' must be an integer.")
             if (burn < 0) stop("'burn' must be non-negative.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("'burn' cannot exceed length of chain.")

             if (length(stat > 1)) stat <- stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("'stat' must be either 'mean' or 'median'")
             if (is.null(stat)) stat <- "mean" # Default to mean

             if (length(method > 1)) method <- method[1]
             if (!(method %in% c("termscore", "prob")))
               stop("'stat' must be either 'termscore' or 'prob'")
             if (is.null(stat)) stat <- "termscore" # Default to termscore

             standardGeneric("gg_coef")
           }
)

#' Generic function to estimate mean/median theta matrix
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("est_theta",
           function(mcmc_fit, burn = 0, thin = 1, stat = "mean") {

             if (!isClass(mcmc_fit, "Lda"))
               stop("'mcmc_fit' must be an Lda object.")
             if (is.null(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if ( (burn %% 1) != 0) stop("'burn' must be an integer.")
             if (burn < 0) stop("'burn' must be non-negative.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("'burn' cannot exceed length of chain.")

             if (length(stat > 1)) stat <- stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("'stat' must be either 'mean' or 'median'")
             if (is.null(stat)) stat <- "mean" # Default to mean

             standardGeneric("est_theta")
           }
)

#' Generic function to estimate mean/median beta matrix
#'
#' @param docs A D x max(N_d) matrix of the words in all documents. Unused word
#'   positions should be set to 0.
#' @export
#' @rdname slda-gettop-methods
setGeneric("est_beta",
           function(mcmc_fit, docs, burn = 0, thin = 1, stat = "mean") {

             if (!isClass(mcmc_fit, "Lda"))
               stop("'mcmc_fit' must be an Lda object.")
             if (is.null(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if ( (burn %% 1) != 0) stop("'burn' must be an integer.")
             if (burn < 0) stop("'burn' must be non-negative.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("'burn' cannot exceed length of chain.")

             if (length(stat > 1)) stat <- stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("'stat' must be either 'mean' or 'median'")
             if (is.null(stat)) stat <- "mean" # Default to mean

             standardGeneric("est_beta")
           }
)
