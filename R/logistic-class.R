#' @include model-class.R
NULL

#' An S4 class to represent a logistic regression model.
#'
#' @slot proposal_sd A (p + 1) x 1 vector of proposal scales for
#'   Metropolis-Hastings sampling of eta.
Logistic <- setClass("Logistic",
  contains = "Model",
  slots = c(
    proposal_sd = "numeric"
  ),
  prototype = list(
    proposal_sd = NA_real_
  )
)

#### proposal_sd
setGeneric("proposal_sd", function(x) standardGeneric("proposal_sd"))

setGeneric("proposal_sd<-", function(x, value) standardGeneric("proposal_sd<-"))

setMethod("proposal_sd", "Logistic", function(x) x@proposal_sd)

setMethod("proposal_sd<-", "Logistic", function(x, value) {
  x@proposal_sd <- value
  x
})

#' Helper function (constructor) for Logistic class.
#'
#' @param proposal_sd A (p + 1) x 1 vector of proposal scales for
#'   Metropolis-Hastings sampling of eta.
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
#' TODO
