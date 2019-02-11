#' Generic function to retrieve the top topics for each document
#'
#' Description of get_toptopics here.
#'
#' @param mcmc_fit An Lda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param ntopics The number of topics to retrieve.
#' @param stat The summary statistic to use on the posterior draws. Default: mean.
#'
#' @return A tibble containing \code{doc}, \code{topic}, and \code{prob}.
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_toptopics",
           function(mcmc_fit, burn = 0, thin = 1, stat = "mean") {

             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")
             if (!isClass(mcmc_fit, "Lda"))
               stop("mcmc_fit must be an Lda object.")
             if ((burn %% 1) != 0) stop("burn must be an integer.")
             if (burn < 0) stop("burn must be non-negative.")
             if ((thin %% 1) != 0) stop("thin must be an integer.")
             if (thin < 1) stop("thin must be positive.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("burn cannot exceed length of chain.")
             if (thin >= (m - burn)) stop("thin cannot exceed length of chain less burn.")

             if (length(stat > 1)) stat = stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("stat must be either 'mean' or 'median'")
             if (is.null(stat)) stat = "mean" # Default to mean

             standardGeneric("get_toptopics")
           }
)

#' Generic function to retrieve the top words for each topic
#'
#' Description of get_topwords here.
#'
#' @param mcmc_fit An Lda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param nwords The number of words to retrieve.
#' @param vocab A character vector containing the vocabulary.
#' @param stat The summary statistic to use on the posterior draws. Default: mean.
#' @param method If "termscore", use term scores (similar to tf-idf). If "prob",
#'   use probabilities. Default: "termscore".
#'
#' @return A tibble containing \code{topic}, \code{word}, and \code{prob}
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_topwords",
           function(beta_, nwords, vocab, method = "termscore") {

             if ((nwords %% 1) != 0) stop("nwords must be an integer.")
             if (nwords < 1) stop("nwords must be an integer greater than 0.")
             if (length(vocab) == 0) stop("vocab must contain at least one element.")
             if (nwords > length(unique(vocab)))
               stop("n_words cannot exceed the number of unique terms in vocab.")

             if (length(method > 1)) method = method[1]
             if (!(method %in% c("termscore", "prob")))
               stop("stat must be either 'termscore' or 'prob'")
             if (is.null(method)) method = "termscore" # Default to termscore

             standardGeneric("get_topwords")
           }
)

#' Generic function to retrieve the empirical topic proportions
#'
#' Compute empirical topic proportions (zbar) from \code{@topics} in an Lda
#' object
#'
#' @param mcmc_fit An Lda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).

setGeneric("get_zbar",
           function(mcmc_fit, burn = 0, thin = 1) {

             if (!isClass(mcmc_fit, "Lda"))
               stop("mcmc_fit must be an Lda object.")
             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")

             if ((burn %% 1) != 0) stop("burn must be an integer.")
             if (burn < 0) stop("burn must be non-negative.")
             if ((thin %% 1) != 0) stop("thin must be an integer.")
             if (thin < 1) stop("thin must be positive.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("burn cannot exceed length of chain.")
             if (thin >= (m - burn)) stop("thin cannot exceed length of chain less burn.")

             standardGeneric("get_zbar")
           }
)

#' Generic function to plot the regression coefficients for sLDA
#'
#' @param mcmc_fit An Lda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param nwords The number of words to retrieve.
#' @param vocab A character vector containing the vocabulary.
#' @param stat The summary statistic to use on the posterior draws. Default: mean.
#' @param method If "termscore", use term scores (similar to tf-idf). If "prob",
#'   use probabilities. Default: "termscore".
#'
#' @return A tibble containing \code{topic}, \code{word}, and \code{prob}
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("gg_coef",
           function(mcmc_fit, beta_, nwords, vocab, varnames, burn = 0, thin = 1,
                    method = "termscore", stat = "mean", errorbw = 0.5) {

             if (!isClass(mcmc_fit, "Lda"))
               stop("mcmc_fit must be an Lda object.")
             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")
             if ((nwords %% 1) != 0) stop("nwords must be an integer.")
             if (nwords < 1) stop("nwords must be an integer greater than 0.")
             if (length(vocab) == 0) stop("vocab must contain at least one element.")
             if (nwords > length(unique(vocab)))
               stop("n_words cannot exceed the number of unique terms in vocab.")
             if ((burn %% 1) != 0) stop("burn must be an integer.")
             if (burn < 0) stop("burn must be non-negative.")
             m <- mcmc_fit@nchain
             if (burn >= m) stop("burn cannot exceed length of chain.")

             if (length(stat > 1)) stat = stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("stat must be either 'mean' or 'median'")
             if (is.null(stat)) stat = "mean" # Default to mean

             if (length(method > 1)) method = method[1]
             if (!(method %in% c("termscore", "prob")))
               stop("stat must be either 'termscore' or 'prob'")
             if (is.null(stat)) stat = "termscore" # Default to termscore

             standardGeneric("gg_coef")
           }
)
