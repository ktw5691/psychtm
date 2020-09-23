#' Compute and visualize the most probable topics in a document or
#' the most probable words in a topic.
#'
#' \code{get_toptopics()} creates a \code{tibble} of topics and the top
#' \code{nwords} words per topic sorted by probability or term score.
#' \code{get_zbar()} computes empirical topic proportions from slot
#' \code{@topics} in an object of class \code{\link{Sldax}}.
#' \code{gg_coef()} plots regression coefficients for \code{\link{Sldax}} models.
#' \code{est_theta()} estimates the mean or median theta matrix.
#' \code{est_beta()} estimates the mean or median beta matrix.
#' \code{get_coherence()} computes the coherence metric for each topic (see Mimno, Wallach, Talley, Leenders, & McCallum, 2011).
#' \code{get_exclusivity()} computes the exclusivity metric for each topic (see Roberts, Stweart, & Airoldi, 2013).
#'
#' @name sldax-gettop-methods
NULL

#' @param mcmc_fit An object of class \code{\link{Sldax}}.
#' @param burn The number of draws to discard as a burn-in period (default: 0).
#' @param thin The number of draws to skip as a thinning period (default: 1; i.e., no thinning).
#' @param stat The summary statistic to use on the posterior draws (default: \code{"mean"}).
#'
#' @rdname sldax-gettop-methods
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

#' @rdname sldax-gettop-methods
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
#'   computed using \code{\link{est_beta()}}. Each row sums to 1.
#' @param docs The \eqn{D} x max(\eqn{N_d}) matrix of documents (word indices)
#'   used to fit the \code{\link{Sldax}} model.
#' @param nwords The number of highest-probability words per topic to consider where
#'   \eqn{M \le V} and \eqn{V} is the size of the corpus vocabulary. (default: \code{10})
#'
#' @rdname sldax-gettop-methods
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

#' @param weight The weight (between 0 and 1) to give to exclusivity (near 1) vs. frequency (near 0). (default: \code{0.7})
#'
#' @rdname sldax-gettop-methods
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
#' @rdname sldax-gettop-methods
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

#' @param nwords The number of words to retrieve (default: V).
#' @param vocab A character vector of length V containing the vocabulary.
#' @param method If "termscore", use term scores (similar to tf-idf). If "prob",
#'   use probabilities (default: "termscore").
#'
#' @return A \eqn{K} x \eqn{V} matrix of term-scores (comparable to tf-idf).
#' @rdname sldax-gettop-methods
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

#' @rdname sldax-gettop-methods
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

#' @param errorbw Positive control parameter for the width of the +/- 2
#'   posterior standard error bars (default: 0.5).
#'
#' @rdname sldax-gettop-methods
#' @export
setGeneric("gg_coef",
           function(mcmc_fit, burn = 0L, thin = 1L, stat = "mean",
                    errorbw = 0.5) {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

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

             if ( !is.numeric(errorbw) ) stop("'errorbw' must be numeric.")
             if ( errorbw <= 0) stop("'errorbw' must be positive.")

             standardGeneric("gg_coef")
           }
)
