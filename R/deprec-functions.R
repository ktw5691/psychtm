#' Generic function to plot regression coefficients for `Sldax` objects
#'
#' @param errorbw Positive control parameter for the width of the +/- 2
#'   posterior standard error bars (default: `0.5`).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' library(lda) # Required if using `prep_docs()`
#' data(teacher_rate)  # Synthetic student ratings of instructors
#' docs_vocab <- prep_docs(teacher_rate, "doc")
#' vocab_len <- length(docs_vocab$vocab)
#' m1 <- gibbs_sldax(rating ~ I(grade - 1), m = 2,
#'                   data = teacher_rate,
#'                   docs = docs_vocab$documents,
#'                   V = vocab_len,
#'                   K = 2,
#'                   model = "sldax")
#' gg_coef(m1)
#' }
#' @rdname sldax-summary
#' @export
setGeneric("gg_coef",
           function(mcmc_fit, burn = 0L, thin = 1L, stat = "mean",
                    errorbw = 0.5) {

             .Deprecated()

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if ( !is.non_negative_integer(burn) )
               stop("'burn' must be a non-negative integer.")

             if ( !is.positive_integer(thin) )
               stop("'thin' must be a positive integer.")

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

#' @importFrom rlang .data
#' @rdname sldax-summary
setMethod("gg_coef",
          c(mcmc_fit = "Sldax"),
          function(mcmc_fit, burn, thin, stat, errorbw) {

            .Deprecated()

            if (!requireNamespace("dplyr", quietly = TRUE)) {
              stop("Package \"dplyr\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            m <- nchain(mcmc_fit)
            keep_index <- seq(burn + 1, m, thin)

            if (stat == "mean") {
              eta <- colMeans(eta(mcmc_fit)[keep_index, ])
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
                  sig = dplyr::if_else(
                    .data$lbcl > 0 | .data$ubcl < 0, "Yes", "No"))),
              sig = factor(.data$sig))

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
