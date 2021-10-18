#' Model class
#'
#' An S4 super class to represent a regression-like model.
#' @rdname Model
#' @slot ndocs The number of documents/observations.
#' @slot nchain The number of iterations of the Gibbs sampler.
#' @slot mu0 A (p + 1) x 1 matrix of prior means for eta.
#' @slot sigma0 A (p + 1) x (p + 1) prior covariance matrix for eta.
#' @slot eta_start A (p + 1) x 1 matrix of starting values for eta.
#' @slot eta A nchain x (p + 1) matrix of draws of regression coefficients.
#' @slot loglike A nchain x 1 vector of the log-likelihood (up to an additive
#'   constant).
#' @slot logpost A nchain x 1 vector of the log-posterior (up to an additive
#'   constant).
#' @slot waic WAIC (up to an additive constant) on the deviance scale.
#' @slot se_waic Standard error of the WAIC.
#' @slot p_eff The effective number of parameters.
#' @slot lpd A nchain x ndocs matrix of predictive posterior likelihoods (NOT
#'   log-likelihoods).
#' @slot extra A list of additional model fitting information. Contains
#'   time_elapsed, start_time, end_time, corrected_label_switching, and call.
Model <- setClass("Model",
  slots = c(
    ndocs       = "numeric",
    nchain      = "numeric",
    mu0         = "numeric",
    sigma0      = "matrix",
    eta_start   = "numeric",
    eta         = "matrix",
    loglike     = "numeric",
    logpost     = "numeric",
    waic        = "numeric",
    se_waic     = "numeric",
    p_eff       = "numeric",
    lpd         = "matrix",
    extra       = "list"),
  prototype = list(
    ndocs = NA_real_,
    nchain = NA_real_,
    mu0 = NA_real_,
    sigma0 = matrix(NA_real_),
    eta_start = NA_real_,
    eta = matrix(NA_real_),
    loglike = NA_real_,
    logpost = NA_real_,
    waic = NA_real_,
    se_waic = NA_real_,
    p_eff = NA_real_,
    lpd = matrix(NA_real_),
    extra = list(time_elapsed = NA_real_,
                 start_time = NA_real_,
                 end_time = NA_real_,
                 corrected_label_switching = FALSE,
                 call = NA_character_)
  )
)

# ndocs
#' @rdname Model
setGeneric("ndocs", function(x) standardGeneric("ndocs"))

#' @rdname Model
setMethod("ndocs", "Model", function(x) x@ndocs)

#' @rdname Model
setGeneric("ndocs<-", function(x, value) standardGeneric("ndocs<-"))

#' @rdname Model
setMethod("ndocs<-", "Model", function(x, value) {
  x@ndocs <- value
  x
})

# nchain
#' @rdname Model
setGeneric("nchain", function(x) standardGeneric("nchain"))

#' @rdname Model
setMethod("nchain", "Model", function(x) x@nchain)

#' @rdname Model
setGeneric("nchain<-", function(x, value) standardGeneric("nchain<-"))

#' @rdname Model
setMethod("nchain<-", "Model", function(x, value) {
  x@nchain <- value
  x
})

# mu0
#' @rdname Model
setGeneric("mu0", function(x) standardGeneric("mu0"))

#' @rdname Model
setMethod("mu0", "Model", function(x) x@mu0)

#' @rdname Model
setGeneric("mu0<-", function(x, value) standardGeneric("mu0<-"))

#' @rdname Model
setMethod("mu0<-", "Model", function(x, value) {
  x@mu0 <- value
  x
})

# sigma0
#' @rdname Model
setGeneric("sigma0", function(x) standardGeneric("sigma0"))

#' @rdname Model
setMethod("sigma0", "Model", function(x) x@sigma0)

#' @rdname Model
setGeneric("sigma0<-", function(x, value) standardGeneric("sigma0<-"))

#' @rdname Model
setMethod("sigma0<-", "Model", function(x, value) {
  x@sigma0 <- value
  x
})

# eta_start
#' @rdname Model
setGeneric("eta_start", function(x) standardGeneric("eta_start"))

#' @rdname Model
setMethod("eta_start", "Model", function(x) x@eta_start)

#' @rdname Model
setGeneric("eta_start<-", function(x, value) standardGeneric("eta_start<-"))

#' @rdname Model
setMethod("eta_start<-", "Model", function(x, value) {
  x@eta_start <- value
  x
})

# eta
#' @rdname Model
setGeneric("eta", function(x) standardGeneric("eta"))

#' @rdname Model
setMethod("eta", "Model", function(x) x@eta)

#' @rdname Model
setGeneric("eta<-", function(x, value) standardGeneric("eta<-"))

#' @rdname Model
setMethod("eta<-", "Model", function(x, value) {
  x@eta <- value
  x
})

# loglike
#' @rdname Model
setGeneric("loglike", function(x) standardGeneric("loglike"))

#' @rdname Model
setMethod("loglike", "Model", function(x) x@loglike)

#' @rdname Model
setGeneric("loglike<-", function(x, value) standardGeneric("loglike<-"))

#' @rdname Model
setMethod("loglike<-", "Model", function(x, value) {
  x@loglike <- value
  x
})

# logpost
#' @rdname Model
setGeneric("logpost", function(x) standardGeneric("logpost"))

#' @rdname Model
setMethod("logpost", "Model", function(x) x@logpost)

setGeneric("logpost<-", function(x, value) standardGeneric("logpost<-"))

#' @rdname Model
setMethod("logpost<-", "Model", function(x, value) {
  x@logpost <- value
  x
})

# waic
#' @rdname Model
setGeneric("waic", function(x) standardGeneric("waic"))

#' @rdname Model
setMethod("waic", "Model", function(x) x@waic)

#' @rdname Model
setGeneric("waic<-", function(x, value) standardGeneric("waic<-"))

#' @rdname Model
setMethod("waic<-", "Model", function(x, value) {
  x@waic <- value
  x
})

# se_waic
#' @rdname Model
setGeneric("se_waic", function(x) standardGeneric("se_waic"))

#' @rdname Model
setMethod("se_waic", "Model", function(x) x@se_waic)

#' @rdname Model
setGeneric("se_waic<-", function(x, value) standardGeneric("se_waic<-"))

#' @rdname Model
setMethod("se_waic<-", "Model", function(x, value) {
  x@se_waic <- value
  x
})

# p_eff
#' @rdname Model
setGeneric("p_eff", function(x) standardGeneric("p_eff"))

#' @rdname Model
setMethod("p_eff", "Model", function(x) x@p_eff)

#' @rdname Model
setGeneric("p_eff<-", function(x, value) standardGeneric("p_eff<-"))

#' @rdname Model
setMethod("p_eff<-", "Model", function(x, value) {
  x@p_eff <- value
  x
})

# lpd
#' @rdname Model
setGeneric("lpd", function(x) standardGeneric("lpd"))

#' @rdname Model
setMethod("lpd", "Model", function(x) x@lpd)

#' @rdname Model
setGeneric("lpd<-", function(x, value) standardGeneric("lpd<-"))

#' @rdname Model
setMethod("lpd<-", "Model", function(x, value) {
  x@lpd <- value
  x
})

# extra
#' @rdname Model
setGeneric("extra", function(x) standardGeneric("extra"))

#' @rdname Model
setMethod("extra", "Model", function(x) x@extra)

#' @rdname Model
setGeneric("extra<-", function(x, value) standardGeneric("extra<-"))

#' @rdname Model
setMethod("extra<-", "Model", function(x, value) {
  x@extra <- value
  x
})

# Helper function (constructor) for Model class.
#' @rdname Model
Model <- function(ndocs, nchain = 1, mu0 = NaN, sigma0 = NaN,
                  eta_start = NaN, eta = NaN, loglike = NaN, logpost = NaN,
                  waic = NaN, se_waic = NaN, p_eff = NaN, lpd = NaN) {
  ndocs <- as.double(ndocs)
  nchain <- as.double(nchain)
  mu0 <- as.double(mu0)
  sigma0 <- as.matrix(sigma0)
  eta_start <- as.double(mu0)
  eta <- as.matrix(eta)
  loglike <- as.double(loglike)
  logpost <- as.double(logpost)
  waic <- as.double(waic)
  se_waic <- as.double(se_waic)
  p_eff <- as.double(p_eff)
  lpd <- as.matrix(lpd)

  new("Model", ndocs = ndocs, nchain = nchain, mu0 = mu0, sigma0 = sigma0,
      eta_start = eta_start, eta = eta, loglike = loglike,
      logpost = logpost, waic = waic, se_waic = se_waic, p_eff = p_eff,
      lpd = lpd)
}

# Validator function for Model class
setValidity("Model", function(object) {
  if ( (object@nchain != NROW(object@eta)) |
       (object@nchain != length(object@loglike)) |
       (object@nchain != length(object@logpost)) |
       (object@nchain != NROW(object@lpd))) {
    "@eta and @lpd must have the same number of rows as the length of @nchain"
  } else if ( (length(object@waic) != length(object@se_waic)) |
              (length(object@waic) != length(object@p_eff))) {
    "@waic, @se_waic, and @p_eff should all be of length 1"
  } else if ( (ncol(object@lpd) != object@ndocs)) {
    "@lpd should have number of columns equal to value of @ndocs"
  } else if ( (length(object@mu0) != length(object@eta_start)) |
              (length(object@mu0) != ncol(object@eta)) |
              (length(object@mu0) != nrow(object@sigma0)) |
              (length(object@mu0) != ncol(object@sigma0))) {
    "@mu0, @eta_start, and @eta should all have the same length which should be equal to the number of rows and columns in @sigma0"
  } else {
    TRUE
  }
})
