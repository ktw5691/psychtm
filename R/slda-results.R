#' @include slda-class.R slda-generic-functions.R
NULL

#' Compute term-scores for each word-topic pair (p. 75, Srivasta & Sahami, 2009)
#'
#' @param beta_ A K x V matrix of V vocabulary probabilities for each of K topics.
#' @return A K x V matrix of term-scores (comparable to tf-idf).
#' @export
term_score = function(beta_) {
  tscore = matrix(0, nrow = nrow(beta_), ncol(beta_))
  for (v in 1:ncol(beta_)) {
    tscore[, v] = beta_[, v] * log(
      exp(log(beta_[, v]) - sum(log(beta_[, v])) / nrow(beta_)))
  }

  return(tscore)
}

#' Return most probable topics for each document.
#'
#' @return A tibble with three columns: \code{doc}: the document number,
#'   \code{topic}: the topic number, \code{prob}: the probability of each topic.
#' @export
setMethod("get_toptopics",
          c(mcmc_fit = "Lda"),
          function(mcmc_fit, burn, thin, stat) {

            m <- mcmc_fit@nchain

            keep_index <- seq(burn + 1, m, thin)
            theta_keep <- mcmc_fit@theta[, , keep_index]

            if (stat == "mean") theta_out <- apply(theta_keep, c(1, 2), mean)
            if (stat == "median") theta_out <- apply(theta_keep, c(1, 2), median)

            ndocs <- mcmc_fit@ndocs
            ntopics <- mcmc_fit@ntopics

            doc_toptopics = matrix(0, ndocs * ntopics, 3)

            row <- 1
            for (d in seq_len(ndocs)) {
              sorted <- sort(theta_out[d, ], decreasing = TRUE)
              for (k in seq_len(ntopics)) {
                doc_toptopics[row + k - 1, 1] <- d
                doc_toptopics[row + k - 1, 2] <- which(
                  theta_out[d, ] == sorted[k])
                doc_toptopics[row + k - 1, 3] <- sorted[k]
              }
              row <- d * ntopics + 1
            }

            doc_toptopics <- tibble::as_tibble(doc_toptopics)
            names(doc_toptopics) <- c("doc", "topic", "prob")
            doc_toptopics$prob <- as.numeric(doc_toptopics$prob)

            return(doc_toptopics)
          }
)

#' Return most probable words for each topic.
#'
#' @param mcmc_fit An Lda object.
#' @return A tibble with three columns: \code{topic}: the topic number,
#'   \code{word}: the vocabulary term, \code{prob}: the term-score or
#'   probability of each word for a given topic.
#' @export
setMethod("get_topwords",
          c(mcmc_fit = "Lda",
            method = "character"),
          function(mcmc_fit, nwords, vocab, burn, thin, method, stat) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)
            beta_keep <- mcmc_fit@beta[, , keep_index]

            if (stat == "mean") beta_out <- apply(beta_keep, c(1, 2), mean)
            if (stat == "median") beta_out <- apply(beta_keep, c(1, 2), median)

            if (method == "termscore") beta_out = term_score(beta_out)

            ntopics <- mcmc_fit@ntopics
            topic_topwords <- matrix("", ntopics * nwords, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(nwords)) {
                topic_topwords[row + v - 1, 1] <- k
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])]
                topic_topwords[row + v - 1, 3] <- sorted[v]
              }
              row <- k * nwords + 1
            }

            topic_topwords <- tibble::as_tibble(topic_topwords)
            names(topic_topwords) <- c("topic", "word", "prob")
            topic_topwords$prob <- as.numeric(topic_topwords$prob)

            return(topic_topwords)
          }
)

#' Return empirical topic proportions
#'
#' Compute empirical topic proportions (zbar) from \code{@topics} in an Lda
#' object
#'
#' @param mcmc_fit An Lda object.
#' @param burn The burn-in period of draws to discard (default: 0).
#' @param thin The thinning period of draws to discard (default: 1; no thinning).
#' @return A D x K matrix of empirical topic proportions (i.e., the relative
#'   frequency of draws for each topic in each document).
#' @export
setMethod("get_zbar",
          c(mcmc_fit = "Lda"),
          function(mcmc_fit, burn, thin) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)

            topics <- mcmc_fit@topics[, , keep_index]
            topics[topics == 0] <- NA # C++ returns 0 if not assigned (no word there)

            # Median topic draw for each word and doc
            z_med <- apply(topics, c(1, 2),
                          function(x) round(median(x, na.rm = TRUE), 0))
            ndoc <- mcmc_fit@ndocs
            ntopic <- mcmc_fit@ntopics

            zbar <- matrix(nrow = ndoc, ncol = ntopic)
            for (d in seq_len(ndoc)) {
              prow <- numeric(ntopic)
              zrow <- z_med[d, ]
              zrow <- zrow[!is.na(zrow)]
              for (j in seq_len(ntopic)) {
                prow[j] <- sum(zrow == j)
              }
              prow <- prow / length(zrow)
              zbar[d, ] <- prow
            }

            return(zbar)
          }
)
