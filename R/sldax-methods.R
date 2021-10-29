#' @include aaa-classes.R aaa-generics.R
NULL

#' @rdname sldax-summary
setMethod("est_beta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat) {

            m <- nchain(mcmc_fit)
            keep   <- seq(burn + 1, m, thin)
            beta_ <- beta_(mcmc_fit)[, , keep]

            if (stat == "mean") beta_hat <- apply(beta_, c(1, 2), mean)
            if (stat == "median") beta_hat <- apply(beta_, c(1, 2), median)

            return(beta_hat)
          }
)

#' @rdname sldax-summary
setMethod("est_theta",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat) {

            m <- nchain(mcmc_fit)
            keep <- seq(burn + 1, m, thin)
            theta <- theta(mcmc_fit)[, , keep]

            if (stat == "mean") theta_hat <- apply(theta, c(1, 2), mean)
            if (stat == "median") theta_hat <- apply(theta, c(1, 2), median)

            return(theta_hat)
          }
)

#' @rdname sldax-summary
setMethod("get_coherence",
          c(beta_ = "matrix", docs = "matrix"),
          function(beta_, docs, nwords) {

            K <- nrow(beta_)
            coher_score <- numeric(K)

            # Get M most probable words for each topic
            topM <- suppressMessages(
              get_topwords(beta_, nwords = nwords, method = "prob"))

            for (k in seq_len(K)) {
              # M most probable words for topic k
              wordsMk_vec <- unlist(topM[topM$topic == k, "word"])

              # Initialize document frequency of word presence
              #   to 0 (max is log(D))
              log_corpusM <- numeric(nwords)
              # Compute log of frequency of topic k's M most frequent words
              #   in corpus
              for (m in seq_len(nwords)) {
                log_corpusM[m] <- log(
                  sum(
                    apply(docs, 1, function(x) wordsMk_vec[m] %in% x)
                  )
                )
              }
              temp_score <- 0
              # Compute coherence for topic k
              for (l in 2:nwords) {
                for (j in seq_len(l - 1)) {
                  # Count presence of co-occurences of two words (1 or 0 for
                  #   each document) in corpus; range is 0 to D
                  d_lj <- sum(
                    apply(
                      docs, 1, function(x) as.integer(
                        (wordsMk_vec[l] %in% x) & (wordsMk_vec[j] %in% x))
                    )
                  )
                  temp_score <- temp_score + log(d_lj + 1) - log_corpusM[j]
                }
              }
              coher_score[k] <- temp_score
            }
          return(coher_score)
          }
)

#' @rdname sldax-summary
setMethod("get_exclusivity",
          c(beta_ = "matrix"),
          function(beta_, nwords, weight) {

            K <- nrow(beta_) # Number of topics
            V <- ncol(beta_) # Vocabulary size

            # Unnormalized marginal word probabilities across topics
            sum_over_topics <- colSums(beta_)
            # Divide each word-topic probability by the marginal word
            #   probabilities (unnormalized)
            ex_temp <- beta_ %*% diag(1 / sum_over_topics) # K x V
            # Exclusivity
            #   V x K matrix of empirical CDF of weighted word probabilities for
            #   each topic
            exc <- apply(ex_temp, 1, rank) / V
            exc <- 1 / exc
            # Frequency
            #   V x K matrix of empirical CDF of word probabilities for each
            #   topic
            fre <- apply(beta_, 1, rank) / V
            fre <- 1 / fre
            # Compute FREX as weighted harmonic mean of frequency and
            #   exclusivity
            #   V x K matrix of FREX values
            frex <- weight * exc + (1 - weight) * fre
            frex <- 1 / frex

            # Get top M word indices per topic
            # M x K matrix
            ind <- apply(beta_, 1, order, decreasing = TRUE)[seq_len(nwords), ]
            # For M most probable words per topic, compute sum of M
            #   corresponding FREX scores
            frex_scores <- numeric(K)
            for (k in seq_len(K)) frex_scores[k] <- sum(frex[ind[, k], k])
            return(frex_scores)
          }
)

#' @rdname sldax-summary
setMethod("get_toptopics",
          c(theta = "matrix"),
          function(theta, ntopics) {

            ndocs   <- NROW(theta)
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

            doc_toptopics <- tibble::as_tibble(doc_toptopics,
                                               .name_repair = "minimal")
            names(doc_toptopics) <- c("doc", "topic", "prob")
            doc_toptopics$prob <- as.numeric(doc_toptopics$prob)

            return(doc_toptopics)
          }
)

#' @rdname sldax-summary
setMethod("get_topwords",
          c(beta_ = "matrix", nwords = "numeric", vocab = "character"),
          function(beta_, nwords, vocab, method) {

            if (method == "termscore") {
              beta_out <- term_score(beta_)
            } else {
              beta_out <- beta_
            }

            ntopics <- NROW(beta_)
            topic_topwords <- matrix("", ntopics * nwords, 3)

            row <- 1
            for (k in seq_len(ntopics)) {
              sorted <- sort(beta_out[k, ], decreasing = TRUE)
              for (v in seq_len(nwords)) {
                topic_topwords[row + v - 1, 1] <- k
                # If tie, pick first match
                topic_topwords[row + v - 1, 2] <- vocab[
                  which(beta_out[k, ] == sorted[v])[1]
                ]
                topic_topwords[row + v - 1, 3] <- sorted[v]
              }
              row <- k * nwords + 1
            }

            topic_topwords <- tibble::as_tibble(topic_topwords,
                                                .name_repair = "minimal")
            names(topic_topwords) <- c("topic", "word", method)
            if (method == "prob") {
              topic_topwords$prob <- as.numeric(topic_topwords$prob)
            } else {
              topic_topwords$termscore <- as.numeric(topic_topwords$termscore)
            }

            return(topic_topwords)
          }
)

#' @rdname sldax-summary
setMethod("get_zbar",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin) {

            m <- nchain(mcmc_fit)
            keep_index <- seq(burn + 1L, m, thin)

            topics <- topics(mcmc_fit)
            if (length(dim(topics)) == 3L) topics <- topics[, , keep_index]

            # C++ returns 0 if not assigned (unused word position)
            topics[topics == 0] <- NA

            # Median topic draw for each word and doc
            if (!is.null(dim(topics)) &&
                length(dim(topics)) == 3L &&
                dim(topics)[3] > 1L) {
              z_med <- apply(
                topics, c(1, 2), function(x) round(median(x, na.rm = TRUE), 0))
            } else {
              z_med <- topics
            }

            ntopic <- ntopics(mcmc_fit)
            if (!is.null(dim(topics)) &&
                length(dim(topics)) == 3 &&
                dim(topics)[3] > 1L) {
              doc_lengths <- apply(topics[, , 1], 1, function(x) sum(!is.na(x)))
              zbar <- t(apply(z_med, 1, tabulate, nbins = ntopic)) / doc_lengths
            } else {
              doc_lengths <- sum(!is.na(topics))
              zbar <- tabulate(z_med, nbins = ntopic) / doc_lengths
            }
            return(zbar)
          }
)

#' @rdname sldax-summary
setMethod("post_regression",
          c(mcmc_fit = "Mlr"),
          function(mcmc_fit) {

            m <- nchain(mcmc_fit)
            burn <- extra(mcmc_fit)$call[["burn"]] # Avoid partial matching
            if (is.null(burn)) burn <- 0L # Default in `gibbs_sldax()` is `0`
            thin <- extra(mcmc_fit)$call[["thin"]] # Avoid partial matching
            if (is.null(thin)) thin <- 1L  # Default in `gibbs_sldax()` is `1`

            mcmc_out <- coda::mcmc(cbind(eta(mcmc_fit), sigma2(mcmc_fit)),
                                   start = burn + 1L, thin = thin)
            colnames(mcmc_out)[ncol(mcmc_out)] <- "sigma2"

            return(mcmc_out)
          }
)

#' @rdname sldax-summary
setMethod("post_regression",
          c(mcmc_fit = "Logistic"),
          function(mcmc_fit) {

            m <- nchain(mcmc_fit)
            burn <- extra(mcmc_fit)$call[["burn"]] # Avoid partial matching
            if (is.null(burn)) burn <- 0L # Default in `gibbs_sldax()` is `0`
            thin <- extra(mcmc_fit)$call[["thin"]] # Avoid partial matching
            if (is.null(thin)) thin <- 1L  # Default in `gibbs_sldax()` is `1`

            return(coda::mcmc(eta(mcmc_fit), start = burn + 1, thin = thin))
          }
)

#' @rdname sldax-summary
setMethod("post_regression",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit) {

            m <- nchain(mcmc_fit)
            burn <- extra(mcmc_fit)$call[["burn"]] # Avoid partial matching
            if (is.null(burn)) burn <- 0L # Default in `gibbs_sldax()` is `0`
            thin <- extra(mcmc_fit)$call[["thin"]] # Avoid partial matching
            if (is.null(thin)) thin <- 1L  # Default in `gibbs_sldax()` is `1`

            # Obtain topic coefficient contrasts
            k <- ntopics(mcmc_fit)
            cm <- ( diag(k, nrow = k) - matrix(rep(1, k ^ 2), nrow = k) ) /
              (k - 1)
            # Number of predictors NOT INCLUDING topics
            n_covs <- ncol(eta(mcmc_fit)) - k
            # Topic "effects" or contrasts
            effectm <- crossprod(
              t(eta(mcmc_fit))[seq(n_covs + 1, n_covs + k), ], cm)
            colnames(effectm) <- paste0("effect_t", seq_len(k))

            if (extra(mcmc_fit)$call$model %in% c("slda", "sldax")) {
              mcmc_out <- coda::mcmc(
                cbind(eta(mcmc_fit), effectm, sigma2(mcmc_fit)),
                start = burn + 1, thin = thin)
              colnames(mcmc_out)[ncol(mcmc_out)] <- "sigma2"
            } else {
              mcmc_out <- coda::mcmc(
                cbind(eta(mcmc_fit), effectm),
                start = burn + 1, thin = thin)
            }

            return(mcmc_out)
          }
)
