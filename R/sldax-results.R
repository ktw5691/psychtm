#' @include sldax-class.R sldax-generic-functions.R
NULL

#' Compute term-scores for each word-topic pair
#'
#' For more details, see Blei, D. M., & Lafferty, J. D. (2009). Topic models. In
#' A. N. Srivastava & M. Sahami (Eds.), Text mining: Classification, clustering,
#' and applications. Chapman and Hall/CRC.
#'
#' @param beta_ A \eqn{K} x \eqn{V} matrix of \eqn{V} vocabulary probabilities
#'   for each of \eqn{K} topics.
#'
#' @return A \eqn{K} x \eqn{V} matrix of term-scores (comparable to tf-idf).
#' @export
term_score <- function(beta_) {

  passed_args <- names(as.list(match.call())[-1])

  if (!"beta_" %in% passed_args)
    stop("Please supply a matrix to 'beta_'.")
  if (!is.matrix(beta_))
    stop("Please supply a matrix to 'beta_'.")
  if (sum(beta_ < 0.0 | beta_ > 1.0) > 0)
    stop("Entries of 'beta_' must be between 0.0 and 1.0.")
  sum_rowsum_beta <- sum(rowSums(beta_))
  K <- nrow(beta_)
  if (sum_rowsum_beta > K + 0.001 | sum_rowsum_beta < K - 0.001)
    stop("Rows of 'beta_' must each sum to 1.0.")

  ntopic <- nrow(beta_)
  ldenom <- apply(log(beta_), 2, sum) / ntopic # Sum logs over topics (rows)
  mdenom <- matrix(ldenom, nrow = ntopic, ncol = ncol(beta_), byrow = TRUE)
  tscore <- beta_ * (log(beta_) - mdenom)

  return(tscore)
}

#' @rdname sldax-gettop-methods
setMethod("get_toptopics",
          c(theta = "matrix"),
          function(theta, ntopics) {

            ndocs   <- nrow(theta)

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

#' @rdname sldax-gettop-methods
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

#' @rdname sldax-gettop-methods
setMethod("get_zbar",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin) {

            m <- nchain(mcmc_fit)
            keep_index <- seq(burn + 1, m, thin)

            topics <- topics(mcmc_fit)[, , keep_index]
            # C++ returns 0 if not assigned (no word there)
            topics[topics == 0] <- NA

            # Median topic draw for each word and doc
            z_med <- apply(topics, c(1, 2),
                           function(x) round(median(x, na.rm = TRUE), 0))

            ntopic <- ntopics(mcmc_fit)
            doc_lengths <- apply(topics[, , 1], 1, function(x) sum(!is.na(x)))

            zbar <- t(apply(z_med, 1, tabulate, nbins = ntopic)) / doc_lengths

            return(zbar)
          }
)

#' @rdname sldax-gettop-methods
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

            m <- nchain(mcmc_fit)
            keep_index <- seq(burn + 1, m, thin)

            if (stat == "mean") {
              eta <- apply(eta(mcmc_fit)[keep_index, ], 2, mean)
            }
            if (stat == "median") {
              eta <- apply(eta(mcmc_fit)[keep_index, ], 2, median)
            }

            varnames <- colnames(eta(mcmc_fit))

            # Regression coefficients for each topic
            coefs <- tibble::tibble(
              coef = eta,
              lbcl = apply(eta(mcmc_fit)[keep_index, ], 2, quantile, .025),
              ubcl = apply(eta(mcmc_fit)[keep_index, ], 2, quantile, .975))

            names(coefs) <- c("est", "lbcl", "ubcl")
            coefs <- cbind(
              coefs, varnames = factor(varnames,
                                       varnames[order(coefs$est, coefs$lbcl)]))
            coefs <- coefs[order(coefs$est), ]

            coefs <- dplyr::mutate(
              dplyr::ungroup(
                dplyr::mutate(
                  dplyr::rowwise(coefs),
                  sig = dplyr::if_else(lbcl > 0 | ubcl < 0, "Yes", "No"))),
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
            return(ggp)
          }
)

#' @rdname sldax-gettop-methods
setMethod("est_theta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat, correct_label_switch, verbose) {

            m <- nchain(mcmc_fit)
            K <- ntopics(mcmc_fit)
            ndoc <- ndocs(mcmc_fit)
            alpha_ <- alpha(mcmc_fit)
            keep <- seq(burn + 1, m, thin)
            topics <- topics(mcmc_fit)[, , keep]
            len <- dim(topics)[3]
            topics[topics == 0] <- NA
            theta <- array(dim = c(ndoc, K, len))

            for (i in seq_len(len)) {
              for (d in seq_len(ndoc)) {
                z_count <- numeric(K)
                for (k in seq_len(K))
                  z_count[k] <- sum(topics[d, , i] == k, na.rm = TRUE)
                theta[d, , i] <- .est_thetad(z_count, alpha_)
              }
              if (verbose) {
                if (!is.nan(i %% floor(len / 10)) & (i %% floor(len / 10) == 0))
                  cat("Iteration", i, "of", len, "\n")
              }
            }

            if (correct_label_switch) {
              # Permute to dimensions: MCMC draws, Docs, Topics
              theta_perm <- aperm(theta, c(3, 1, 2))
              # Relabeling algorithm from Stephens (2000)
              relabel_out <- label.switching::stephens(theta_perm)
              if (verbose)
                cat("Relabeling algorithm status: ", relabel_out$status, "\n")
              reorder_theta <- label.switching::permute.mcmc(
                aperm(theta_perm, c(1, 3, 2)),
                permutations = relabel_out$permutations)$output
              reorder_theta <- aperm(reorder_theta, c(1, 3, 2)) # MCMC, D, K
              if (stat == "mean")
                theta_mean <- apply(reorder_theta, c(2, 3), mean)
              if (stat == "median")
                theta_mean <- apply(reorder_theta, c(2, 3), median)
            } else {
              if (stat == "mean") theta_mean <- apply(theta, c(1, 2), mean)
              if (stat == "median") theta_mean <- apply(theta, c(1, 2), median)
            }
            return(theta_mean)
          }
)

#' @rdname sldax-gettop-methods
setMethod("est_beta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, docs, burn, thin, stat,
                   correct_label_switch, verbose) {

            m <- nchain(mcmc_fit)
            K <- ntopics(mcmc_fit)
            V <- nvocab(mcmc_fit)
            gamma_ <- gamma(mcmc_fit)
            keep   <- seq(burn + 1, m, thin)
            topics <- topics(mcmc_fit)[, , keep]
            len    <- dim(topics)[3]

            topics[topics == 0] <- NA
            beta_ <- array(dim = c(K, V, len))

            for (i in seq_len(len)) {
              wz_co <- .count_topic_word(K, V, topics[, , i], docs)
              for (k in seq_len(K)) beta_[k, , i] <- .est_betak(
                wz_co[k, ], gamma_)
              if (verbose) {
                if (!is.nan(i %% floor(len / 10)) & (i %% floor(len / 10) == 0))
                  cat("Iteration", i, "of", len, "\n")
              }
            }

            if (correct_label_switch) {
              # Permute to dimensions: MCMC draws, Words, Topics
              beta_perm <- aperm(beta_, c(3, 2, 1))
              # Relabeling algorithm from Stephens (2000)
              relabel_out <- label.switching::stephens(beta_perm)
              if (verbose)
                cat("Relabeling algorithm status: ", relabel_out$status, "\n")
              reorder_beta <- label.switching::permute.mcmc(
                aperm(beta_perm, c(1, 3, 2)),
                permutations = relabel_out$permutations)$output
              reorder_beta <- aperm(reorder_beta, c(1, 3, 2)) # MCMC, K, V
              if (stat == "mean")
                beta_mean <- t(apply(reorder_beta, c(2, 3), mean))
              if (stat == "median")
                beta_mean <- t(apply(reorder_beta, c(2, 3), median))
            } else {
              if (stat == "mean") beta_mean <- apply(beta_, c(1, 2), mean)
              if (stat == "median") beta_mean <- apply(beta_, c(1, 2), median)
            }
            return(beta_mean)
          }
)
