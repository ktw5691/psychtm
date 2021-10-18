#' @include model-class.R
NULL

#' Logistic class
#'
#' S4 class for a logistic regression model. Inherits from \code{\linkS4class{Model}}.
#' @name Logistic
#' @rdname Logistic
#' @slot proposal_sd A vector of p + 1 proposal scales/standard deviations for
#'   sampling of p + 1 regression coefficients by Metropolis-Hastings.
Logistic <- setClass("Logistic",
  contains = "Model",
  slots = c(
    proposal_sd = "numeric"
  ),
  prototype = list(
    proposal_sd = NA_real_
  )
)

# proposal_sd slot generic
#' @rdname Logistic
setGeneric("proposal_sd", function(x) standardGeneric("proposal_sd"))

# Slot \code{@proposal_sd} generic accessor
#' @rdname Logistic
setMethod("proposal_sd", "Logistic", function(x) x@proposal_sd)

# Slot \code{@proposal_sd} generic setter function
#' @rdname Logistic
setGeneric("proposal_sd<-", function(x, value) standardGeneric("proposal_sd<-"))

# Slot \code{@proposal_sd} setter method
#' @rdname Logistic
setMethod("proposal_sd<-", "Logistic", function(x, value) {
  x@proposal_sd <- value
  validObject(x)
  x
})

# Logistic helper function (constructor) for Logistic class.
#' @rdname Logistic
Logistic <- function(proposal_sd = NaN, ...) {
  super <- Model(...)
  proposal_sd <- as.double(proposal_sd)

  new("Logistic", proposal_sd = proposal_sd,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

#' Validator function for Logistic class
#' @rdname Logistic
#' TODO
