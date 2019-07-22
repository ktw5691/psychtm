#' @include slda-class.R slda-generic-functions.R
NULL

#' Compute term-scores for each word-topic pair (p. 75, Srivasta & Sahami, 2009)
#'
#' @param beta_ A \eqn{K} x \eqn{V} matrix of \eqn{V} vocabulary probabilities
#'   for each of \eqn{K} topics.
#'
#' @return A \eqn{K} x \eqn{V} matrix of term-scores (comparable to tf-idf).
#' @export
term_score <- function(beta_) {

  tscore <- matrix(0, nrow = nrow(beta_), ncol(beta_))
  for (v in 1:ncol(beta_)) {
    tscore[, v] <- beta_[, v] * log(
      exp(log(beta_[, v]) - sum(log(beta_[, v])) / nrow(beta_)))
  }

  return(tscore)
}

setMethod("get_toptopics",
          c(theta = "matrix", ntopics = "numeric"),
          function(theta, ntopics) {

            ndocs   <- nrow(theta)
            ntopics <- ncol(theta)

            doc_toptopics <- matrix(0, ndocs * ntopics, 3)

            row <- 1
            for (d in seq_len(ndocs)) {
              sorted <- sort(theta[d, ], decreasing = TRUE)
              for (k in seq_len(ntopics)) {
                doc_toptopics[row + k - 1, 1] <- d
                doc_toptopics[row + k - 1, 2] <- which(
                  theta[d, ] == sorted[k])[1] # If tie, pick first match
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

setMethod("get_topwords",
          c(beta_ = "matrix", nwords = "numeric", vocab = "character"),
          function(beta_, nwords, vocab, method) {

            if (method == "termscore") {
              beta_out <- term_score(beta_)
            } else {
              beta_out <- beta_
            }

            ntopics <- nrow(beta_)
            topic_topwords <- matrix("", ntopics * nwords, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(nwords)) {
                topic_topwords[row + v - 1, 1] <- k
                # If tie, pick first match
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])[1]]
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

setMethod("get_zbar",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin) {

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)

            topics <- mcmc_fit@topics[, , keep_index]
            # C++ returns 0 if not assigned (no word there)
            topics[topics == 0] <- NA

            # Median topic draw for each word and doc
            z_med <- apply(topics, c(1, 2),
                           function(x) round(median(x, na.rm = TRUE), 0))

            ntopic <- mcmc_fit@ntopics
            doc_lengths <- apply(topics[, , 1], 1, function(x) sum(!is.na(x)))

            zbar <- t(apply(z_med, 1, tabulate, nbins = ntopic)) / doc_lengths

            return(zbar)
          }
)

setMethod("gg_coef",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat, errorbw) {

            if (!requireNamespace("dplyr", quietly = TRUE)) {
              stop("Package \"dplyr\" needed for this function to work.
                    Please install it.",
                   call. = FALSE)
            }
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work.
                    Please install it.",
                   call. = FALSE)
            }

            m <- mcmc_fit@nchain
            keep_index <- seq(burn + 1, m, thin)

            if (stat == "mean") {
              eta <- apply(mcmc_fit@eta[keep_index, ], 2, mean)
            }
            if (stat == "median") {
              eta <- apply(mcmc_fit@eta[keep_index, ], 2, median)
            }

            varnames <- colnames(mcmc_fit@eta)

            # Regression coefficients for each topic
            coefs <- tibble::tibble(
              coef = eta,
              lbcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .025),
              ubcl = apply(mcmc_fit@eta[keep_index, ], 2, quantile, .975))

            names(coefs) <- c("est", "lbcl", "ubcl")
            coefs <- cbind(
              coefs, varnames = factor(varnames,
                                       varnames[order(coefs$est, coefs$lbcl)]))
            coefs <- coefs[order(coefs$est), ]

            coefs <- dplyr::mutate(
              dplyr::ungroup(
                dplyr::mutate(
                  dplyr::rowwise(coefs),
                  sig = dplyr::if_else(lbcl > 0 || ubcl < 0, "Yes", "No"))),
              sig = factor(sig))

            ggp <- ggplot2::ggplot(
              coefs, ggplot2::aes_string(x = "varnames", y = "est",
                                         color = "sig"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_point())
            ggp <- ggplot2::`%+%`(
              ggp, ggplot2::geom_errorbar(
                width = errorbw, ggplot2::aes_string(ymin = "lbcl",
                                                     ymax = "ubcl")))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::geom_hline(yintercept = 0))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::theme_bw())
            ggp <- ggplot2::`%+%`(
              ggp, ggplot2::scale_color_discrete("Non-Zero"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::xlab("Predictor"))
            ggp <- ggplot2::`%+%`(ggp, ggplot2::ylab("Posterior Estimate"))
            ggp <- ggplot2::`%+%`(
              ggp, ggplot2::theme(axis.text.x = ggplot2::element_text(
                size = 10, angle = 57.5, hjust = 1), legend.position = "left"))
            print(ggp)
          }
)

setMethod("est_theta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat) {

            m <- mcmc_fit@nchain
            K <- mcmc_fit@ntopics
            ndoc <- mcmc_fit@ndocs
            alpha_ <- mcmc_fit@alpha
            keep <- seq(burn + 1, m, thin)
            topics <- mcmc_fit@topics[, , keep]
            len <- dim(topics)[3]
            topics[topics == 0] <- NA
            theta <- array(dim = c(ndoc, ncol = K, len))

            for (i in seq_len(len)) {
              for (d in seq_len(ndoc)) {
                z_count <- numeric(K)
                for (k in seq_len(K))
                  z_count[k] <- sum(topics[d, , i] == k, na.rm = TRUE)
                theta[d, , i] <- .est_thetad(z_count, alpha_)
              }
              if (i %% floor(len / 10) == 0)
                cat("Iteration", i, "of", len, "\n")
            }

            if (stat == "mean") theta_mean <- apply(theta, c(1, 2), mean)
            if (stat == "median") theta_mean <- apply(theta, c(1, 2), median)
            return(theta_mean)
          }
)

setMethod("est_beta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, docs, burn, thin, stat) {

            m <- mcmc_fit@nchain
            K <- mcmc_fit@ntopics
            V <- mcmc_fit@nvocab
            ndoc <- mcmc_fit@ndocs
            gamma_ <- mcmc_fit@gamma
            keep <- seq(burn + 1, m, thin)
            topics <- mcmc_fit@topics[, , keep]
            len <- dim(topics)[3]
            topics[topics == 0] <- NA
            beta_ <- array(dim = c(K, ncol = V, len))

            for (i in seq_len(len)) {
              wz_co <- .count_topic_word(K, V, topics[, , i], docs)
              for (k in seq_len(K)) beta_[k, , i] <- .est_betak(
                wz_co[k, ], gamma_)
              if (i %% floor(len / 10) == 0)
                cat("Iteration", i, "of", len, "\n")
            }

            if (stat == "mean") beta_mean <- apply(beta_, c(1, 2), mean)
            if (stat == "median") beta_mean <- apply(beta_, c(1, 2), median)
            return(beta_mean)
          }
)
