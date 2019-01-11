#' Generic function to retrieve the top topics for each document
#'
#' Description of get_toptopics here.
#'
#' @param mcmc_fit An SLda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param n_topics The number of topics to retrieve.
#' @param stat The summary statistic to use on the posterior draws. Default: mean.
#'
#' @return A tibble containing \code{doc}, \code{topic}, and \code{prob}.
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_toptopics",
           function(mcmc_fit, n_topics, burn = 0, thin = 1, stat = "mean") {

             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")
             if (!isClass(mcmc_fit, "SLda"))
               stop("mcmc_fit must be an SLda object.")
             if ((n_topics %% 1) != 0) stop("n_topics must be an integer.")
             if (n_topics < 1) stop("n_topics must be an integer greater than 0.")
             K_ <- mcmc_fit@ntopics
             if (n_topics > K_) stop(paste0("n_topics cannot exceed the number of topics in the estimated model: ", K_))
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
#' @param mcmc_fit An SLda object.
#' @param burn The number of draws to discard as a burn-in period. Default: 0.
#' @param thin The number of draws to skip as a thinning period. Default: 1 (no thinning).
#' @param n_words The number of words to retrieve.
#' @param vocab A character vector containing the vocabulary.
#' @param stat The summary statistic to use on the posterior draws. Default: mean.
#'
#' @return A tibble containing \code{topic}, \code{word}, and \code{prob}
#'
#' @export
#' @rdname slda-gettop-methods
setGeneric("get_topwords",
           function(mcmc_fit, n_words, vocab, burn = 0, thin = 1,
                    stat = "mean") {

             if (!isClass(mcmc_fit, "SLda"))
               stop("mcmc_fit must be an SLda object.")
             if (is.null(mcmc_fit)) stop("Please supply an object to mcmc_fit.")
             if ((n_words %% 1) != 0) stop("n_words must be an integer.")
             if (n_words < 1) stop("n_words must be an integer greater than 0.")
             if (length(vocab) == 0) stop("vocab must contain at least one element.")
             if (n_words > length(unique(vocab))) stop("n_words cannot exceed the number of unique terms in vocab.")
             if ((burn %% 1) != 0) stop("burn must be an integer.")
             if (burn < 0) stop("burn must be non-negative.")
             if ((thin %% 1) != 0) stop("thin must be an integer.")
             if (thin < 1) stop("thin must be positive.")
             m <-  mcmc_fit@nchain
             if (burn >= m) stop("burn cannot exceed length of chain.")
             if (thin >= (m - burn)) stop("thin cannot exceed length of chain less burn.")

             if (length(stat > 1)) stat = stat[1]
             if (!(stat %in% c("mean", "median")))
               stop("stat must be either 'mean' or 'median'")
             if (is.null(stat)) stat = "mean" # Default to mean

             standardGeneric("get_topwords")
           }
)
