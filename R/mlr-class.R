#' @include model-class.R
NULL

#' Mlr class
#'
#' S4 class for a regression model. Inherits from \code{\linkS4class{Model}}.
#' @rdname Mlr
#' @slot a0 A prior shape hyperparameter for sigma2.
#' @slot b0 A prior rate hyperparameter for sigma2.
#' @slot sigma2 A nchain x 1 numeric vector of draws of the residual variance.
Mlr <- setClass("Mlr",
  contains = "Model",
  slots = c(
    a0          = "numeric",
    b0          = "numeric",
    sigma2      = "numeric"),
  prototype = list(
    a0 = NA_real_,
    b0 = NA_real_,
    sigma2 = NA_real_
  )
)

#' @rdname Mlr
setGeneric("sigma2", function(x) standardGeneric("sigma2"))

#' @rdname Mlr
setMethod("sigma2", "Mlr", function(x) x@sigma2)

#' @rdname Mlr
setGeneric("sigma2<-", function(x, value) standardGeneric("sigma2<-"))

#' @rdname Mlr
setMethod("sigma2<-", "Mlr", function(x, value) {
  x@sigma2 <- value
  x
})

#' @rdname Mlr
setGeneric("a0", function(x) standardGeneric("a0"))

#' @rdname Mlr
setMethod("a0", "Mlr", function(x) x@a0)

#' @rdname Mlr
setGeneric("a0<-", function(x, value) standardGeneric("a0<-"))

#' @rdname Mlr
setMethod("a0<-", "Mlr", function(x, value) {
  x@a0 <- value
  x
})

#' @rdname Mlr
setGeneric("b0", function(x) standardGeneric("b0"))

#' @rdname Mlr
setMethod("b0", "Mlr", function(x) x@b0)

#' @rdname Mlr
setGeneric("b0<-", function(x, value) standardGeneric("b0<-"))

#' @rdname Mlr
setMethod("b0<-", "Mlr", function(x, value) {
  x@b0 <- value
  x
})

# Helper function (constructor) for Mlr class.
#' @rdname Mlr
Mlr <- function(a0 = 0.001, b0 = 0.001, sigma2 = NaN, ...) {
  super <- Model(...)
  a0 <- as.double(a0)
  b0 <- as.double(b0)
  sigma2 <- as.double(sigma2)

  new("Mlr", a0 = a0, b0 = b0, sigma2 = sigma2,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

#' Validator function for Mlr class
#' @rdname Mlr
#' TODO
