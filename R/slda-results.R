#' @include slda-class.R slda-generic-functions.R
NULL

#' Return most probable words for topic
#' @export
setMethod("get_toptopics",
          c(mcmc_fit = "Lda", n_topics = "numeric"),
          function(mcmc_fit, n_topics, burn, thin, stat) {

            m <- mcmc_fit@nchain

            keep_index <- seq(burn + 1, m, thin)
            theta_keep <- mcmc_fit@theta[, , keep_index]

            if (stat == "mean") theta_out <- apply(theta_keep, c(1, 2), mean)
            if (stat == "median") theta_out <- apply(theta_keep, c(1, 2), median)

            ndocs = mcmc_fit@ndocs

            doc_toptopics = matrix(0, ndocs * n_topics, 3)

            row <- 1
            for (d in seq_len(ndocs)) {
              sorted <- sort(theta_out[d, ], decreasing = TRUE)
              for (k in seq_len(n_topics)) {
                doc_toptopics[row + k - 1, 1] <- d
                doc_toptopics[row + k - 1, 2] <- which(
                  theta_out[d, ] == sorted[k])
                doc_toptopics[row + k - 1, 3] <- sorted[k]
              }
              row <- d * n_topics + 1
            }

            doc_toptopics <- tibble::as_tibble(doc_toptopics)
            names(doc_toptopics) <- c("doc", "topic", "prob")
            doc_toptopics$prob <- as.numeric(doc_toptopics$prob)

            return(doc_toptopics)
          }
)

#' Return most probable words for document
#' @export
setMethod("get_topwords",
          c(mcmc_fit = "Lda", n_words = "numeric", vocab = "character"),
          function(mcmc_fit, n_words, vocab, burn, thin, stat) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)
            beta_keep <- mcmc_fit@beta[, , keep_index]

            if (stat == "mean") beta_out <- apply(beta_keep, c(1, 2), mean)
            if (stat == "median") beta_out <- apply(beta_keep, c(1, 2), median)

            ntopics <- mcmc_fit@ntopics
            topic_topwords <- matrix("", ntopics * n_words, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(n_words)) {
                topic_topwords[row + v - 1, 1] <- k
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])]
                topic_topwords[row + v - 1, 3] <- sorted[v]
              }
              row <- k * n_words + 1
            }

            topic_topwords <- tibble::as_tibble(topic_topwords)
            names(topic_topwords) <- c("topic", "word", "prob")
            topic_topwords$prob <- as.numeric(topic_topwords$prob)

            return(topic_topwords)
          }
)

#' Return most probable words for document
#' @export
setMethod("get_topwords",
          c(mcmc_fit = "Lda", n_words = "numeric", vocab = "numeric"),
          function(mcmc_fit, n_words, vocab, burn, thin, stat) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)
            beta_keep <- mcmc_fit@beta[, , keep_index]

            if (stat == "mean") beta_out <- apply(beta_keep, c(1, 2), mean)
            if (stat == "median") beta_out <- apply(beta_keep, c(1, 2), median)

            ntopics <- mcmc_fit@ntopics
            topic_topwords <- matrix("", ntopics * n_words, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(n_words)) {
                topic_topwords[row + v - 1, 1] <- k
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])]
                topic_topwords[row + v - 1, 3] <- sorted[v]
              }
              row <- k * n_words + 1
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
#' @export
setMethod("get_zbar",
          c(mcmc_fit = "Lda"),
          function(mcmc_fit, burn, thin) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)
            beta_keep <- mcmc_fit@beta[, , keep_index]

            # Exclude 0 ("missing value" in C++ code for unused word positions)
            zbarm = apply(mcmc_fit@topics, 3, function(x) table(x, exclude = 0))
            zbarm = apply(zbarm, 2, function(x) x / sum(x))

            return(zbarm)
          }
)
