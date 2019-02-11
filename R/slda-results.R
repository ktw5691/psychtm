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
                  theta_out[d, ] == sorted[k])[1] # If tie, pick first match
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
#' @param beta_ A K x V matrix of topic-word probabilities.
#' @return A tibble with three columns: \code{topic}: the topic number,
#'   \code{word}: the vocabulary term, \code{prob}: the term-score or
#'   probability of each word for a given topic.
setMethod("get_topwords",
          c(beta_ = "matrix", nwords = "numeric", vocab = "character"),
          function(beta_, nwords, vocab, method) {

            if (method == "termscore") {
              beta_out = term_score(beta_)
            } else {
              beta_out = beta_
            }

            ntopics <- nrow(beta_)
            topic_topwords <- matrix("", ntopics * nwords, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(nwords)) {
                topic_topwords[row + v - 1, 1] <- k
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])[1]] # If tie, pick first match
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

setMethod("gg_coef",
          c(mcmc_fit = "Slda"),
          function(mcmc_fit, beta_, nwords, vocab, varnames, burn, thin, method, stat, errorbw) {

            if (!requireNamespace("dplyr", quietly = TRUE)) {
              stop("Package \"dplyr\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)

            if (stat == "mean") {
              eta = apply(mcmc_fit@eta[keep_index, ], 2, mean)
            }
            if (stat == "median") {
              eta = apply(mcmc_fit@eta[keep_index, ], 2, median)
            }

            topic_dist = get_topwords(beta_, nwords, vocab, method)

            len = ncol(mcmc_fit@eta)
            ntopics = mcmc_fit@ntopics
            # Top words per topic
            topics_top = vector("character", len)
            if (is.null(varnames) | length(varnames) != len - ntopics) {
              topics_top[seq_len(len - ntopics)] <- paste0("V", seq_len(len - ntopics))
            } else {
              topics_top[seq_len(len - ntopics)] <- varnames
            }
            for (k in seq_len(ntopics)) {
              temp = dplyr::filter(topic_dist, topic == k)
              temp = dplyr::top_n(temp, nwords, prob)
              temp = paste(temp$word, collapse = ", ")
              topics_top[len - ntopics + k] = temp
            }

            # Regression coefficients for each topic
            coefs <- tibble::tibble(
              coef = eta,
              lbcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .025),
              ubcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .975))

            names(coefs) = c("est", "lbcl", "ubcl")
            coefs <- cbind(
              coefs, topics = factor(topics_top,
                                     topics_top[order(coefs$est, coefs$lbcl)]))
            coefs <- coefs[order(coefs$est), ]

            coefs <- dplyr::mutate(
              dplyr::ungroup(
                dplyr::mutate(dplyr::rowwise(coefs),
                              sig = dplyr::if_else(lbcl > 0 || ubcl < 0, 1, 0))),
              sig = factor(sig))

            ggp <- ggplot2::ggplot(coefs, ggplot2::aes(topics, est, color = sig))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_point())
            ggp <- ggplot2::`%+%`(ggp,
                                  ggplot2::geom_errorbar(width = errorbw,
                                                         ggplot2::aes(ymin = lbcl,
                                                                      ymax = ubcl)))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_hline(yintercept = 0))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::theme_bw())
            ggp <- ggplot2::`%+%`(ggp, ggplot2::scale_color_discrete("Non-Zero"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::xlab("Predictor"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::ylab("Posterior Estimate"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::theme(axis.text.x = ggplot2::element_text(
              size = 10, angle = 57.5, hjust = 1), legend.position = "left"))
            print(ggp)
          }
)

setMethod("gg_coef",
          c(mcmc_fit = "Sldalogit"),
          function(mcmc_fit, beta_, nwords, vocab, varnames, burn, thin, method, stat, errorbw) {

            if (!requireNamespace("dplyr", quietly = TRUE)) {
              stop("Package \"dplyr\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)

            if (stat == "mean") {
              eta = apply(mcmc_fit@eta[keep_index, ], 2, mean)
            }
            if (stat == "median") {
              eta = apply(mcmc_fit@eta[keep_index, ], 2, median)
            }

            topic_dist = get_topwords(beta_, nwords, vocab, method)

            len = ncol(mcmc_fit@eta)
            ntopics = mcmc_fit@ntopics
            # Top words per topic
            topics_top = vector("character", len)
            if (is.null(varnames) | length(varnames) != len - ntopics) {
              topics_top[seq_len(len - ntopics)] <- paste0("V", seq_len(len - ntopics))
            } else {
              topics_top[seq_len(len - ntopics)] <- varnames
            }
            for (k in seq_len(ntopics)) {
              temp = dplyr::filter(topic_dist, topic == k)
              temp = dplyr::top_n(temp, nwords, prob)
              temp = paste(temp$word, collapse = ", ")
              topics_top[len - ntopics + k] = temp
            }

            # Regression coefficients for each topic
            coefs <- tibble::tibble(
              coef = eta,
              lbcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .025),
              ubcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .975))

            names(coefs) = c("est", "lbcl", "ubcl")
            coefs <- cbind(
              coefs, topics = factor(topics_top,
                                     topics_top[order(coefs$est, coefs$lbcl)]))
            coefs <- coefs[order(coefs$est), ]

            coefs <- dplyr::mutate(
              dplyr::ungroup(
                dplyr::mutate(dplyr::rowwise(coefs),
                              sig = dplyr::if_else(lbcl > 0 || ubcl < 0, 1, 0))),
              sig = factor(sig))

            ggp <- ggplot2::ggplot(coefs, ggplot2::aes(topics, est, color = sig))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_point())
            ggp <- ggplot2::`%+%`(ggp,
                                  ggplot2::geom_errorbar(width = errorbw,
                                                         ggplot2::aes(ymin = lbcl,
                                                                      ymax = ubcl)))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_hline(yintercept = 0))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::theme_bw())
            ggp <- ggplot2::`%+%`(ggp, ggplot2::scale_color_discrete("Non-Zero"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::xlab("Predictor"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::ylab("Posterior Estimate"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::theme(axis.text.x = ggplot2::element_text(
              size = 10, angle = 57.5, hjust = 1), legend.position = "left"))
            print(ggp)
          }
)
