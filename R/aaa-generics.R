#' @include aaa-classes.R
NULL

########----------------------------------------------------------------########
########                            Model Class                         ########
########----------------------------------------------------------------########

#' Create generic `ndocs` function for class
#'
#' @param x An `Model` object.
#'
#' @return Integer number of documents.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' ndocs(m1)
#' @export
setGeneric("ndocs", function(x) standardGeneric("ndocs"))

#' @rdname Model-class
setMethod("ndocs", "Model", function(x) x@ndocs)

#' Create generic `ndocs<-` function for class
#'
#' @param x An `Model` object.
#' @param value Integer number of documents to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' ndocs(m1) <- 2
#' @export
setGeneric("ndocs<-", function(x, value) standardGeneric("ndocs<-"))

#' @rdname Model-class
setMethod("ndocs<-", "Model", function(x, value) {
  x@ndocs <- value
  x
})

#' Create generic `nchain` function for class
#'
#' @param x An `Model` object.
#'
#' @return Integer length of sampler chain.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' nchain(m1)
#' @export
setGeneric("nchain", function(x) standardGeneric("nchain"))

#' @rdname Model-class
setMethod("nchain", "Model", function(x) x@nchain)

#' Create generic `nchain<-` function for class
#'
#' @param x An `Model` object.
#' @param value Integer length of sampler chain to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' nchain(m1) <- 100
#' @export
setGeneric("nchain<-", function(x, value) standardGeneric("nchain<-"))

#' @rdname Model-class
setMethod("nchain<-", "Model", function(x, value) {
  x@nchain <- value
  x
})

#' Create generic `mu0` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric vector of prior means for regression coefficients `eta`.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' mu0(m1)
#' @export
setGeneric("mu0", function(x) standardGeneric("mu0"))

#' @rdname Model-class
setMethod("mu0", "Model", function(x) x@mu0)

#' Create generic `mu0<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric vector of prior means for regression coefficients to
#'   assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' mu0(m1) <- rep(0.0, times = 2)
#' @export
setGeneric("mu0<-", function(x, value) standardGeneric("mu0<-"))

#' @rdname Model-class
setMethod("mu0<-", "Model", function(x, value) {
  x@mu0 <- value
  x
})

#' Create generic `sigma0` function for class
#'
#' @param x An `Model` object.
#'
#' @return Double matrix of prior variances and covariances for regression
#'   coefficients.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' sigma0(m1)
#' @export
setGeneric("sigma0", function(x) standardGeneric("sigma0"))

#' @rdname Model-class
setMethod("sigma0", "Model", function(x) x@sigma0)

#' Create generic `sigma0<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric covariance matrix of prior for regression coefficients
#'   to assign to slot.
#'
#' @return None.
#' @examples
#' m1 <- Model(ndocs = 1)
#' sigma0(m1) <- diag(1.0, 2)
#'
#' @export
setGeneric("sigma0<-", function(x, value) standardGeneric("sigma0<-"))

#' @rdname Model-class
setMethod("sigma0<-", "Model", function(x, value) {
  x@sigma0 <- value
  x
})

#' Create generic `eta_start` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric vector of starting values for regression coefficients.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' eta_start(m1)
#' @export
setGeneric("eta_start", function(x) standardGeneric("eta_start"))

#' @rdname Model-class
setMethod("eta_start", "Model", function(x) x@eta_start)

#' Create generic `eta_start<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric vector of starting values for regression coefficients to
#'   assign to slot.
#'
#' @return None.
#' @examples
#' m1 <- Model(ndocs = 1)
#' eta_start(m1) <- rep(0.0, times = 2)
#'
#' @export
setGeneric("eta_start<-", function(x, value) standardGeneric("eta_start<-"))

#' @rdname Model-class
setMethod("eta_start<-", "Model", function(x, value) {
  x@eta_start <- value
  x
})

#' Create generic `eta` function for class
#'
#' @param x An `Model` object.
#'
#' @return A numeric matrix of posterior draws of regression coefficients.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' eta(m1)
#' @export
setGeneric("eta", function(x) standardGeneric("eta"))

#' @rdname Model-class
setMethod("eta", "Model", function(x) x@eta)

#' Create generic `eta<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric vector of regression coefficients to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' eta(m1) <- matrix(c(-1.0, 1.0), nrow = 1, ncol = 2)
#' @export
setGeneric("eta<-", function(x, value) standardGeneric("eta<-"))

#' @rdname Model-class
setMethod("eta<-", "Model", function(x, value) {
  x@eta <- value
  x
})

#' Create generic `loglike` function for class
#'
#' @param x An `Model` object.
#'
#' @return A numeric vector of log-likelihood values across sampler iterations.
#' @examples
#' m1 <- Model(ndocs = 1)
#' loglike(m1)
#'
#' @export
setGeneric("loglike", function(x) standardGeneric("loglike"))

#' @rdname Model-class
setMethod("loglike", "Model", function(x) x@loglike)

#' Create generic `loglike<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric vector of log likelihoods to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' loglike(m1) <- rep(NaN, times = nchain(m1))
#' @export
setGeneric("loglike<-", function(x, value) standardGeneric("loglike<-"))

#' @rdname Model-class
setMethod("loglike<-", "Model", function(x, value) {
  x@loglike <- value
  x
})

#' Create generic `logpost` function for class
#'
#' @param x An `Model` object.
#'
#' @return A numeric vector of log-posterior values across sampler iterations.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' logpost(m1)
#' @export
setGeneric("logpost", function(x) standardGeneric("logpost"))

#' @rdname Model-class
setMethod("logpost", "Model", function(x) x@logpost)

#' Create generic `logpost<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric vector of log posteriors to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' logpost(m1) <- rep(NaN, times = nchain(m1))
#' @export
setGeneric("logpost<-", function(x, value) standardGeneric("logpost<-"))

#' @rdname Model-class
setMethod("logpost<-", "Model", function(x, value) {
  x@logpost <- value
  x
})

#' Create generic `waic` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric value of the Watanabe Information Criterion (WAIC).
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' waic(m1)
#' @export
setGeneric("waic", function(x) standardGeneric("waic"))

#' @rdname Model-class
setMethod("waic", "Model", function(x) x@waic)

#' Create generic `waic<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric value of WAIC to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' waic(m1) <- NaN
#' @export
setGeneric("waic<-", function(x, value) standardGeneric("waic<-"))

#' @rdname Model-class
setMethod("waic<-", "Model", function(x, value) {
  x@waic <- value
  x
})

#' Create generic `se_waic` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric standard error for WAIC estimate.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' se_waic(m1)
#' @export
setGeneric("se_waic", function(x) standardGeneric("se_waic"))

#' @rdname Model-class
setMethod("se_waic", "Model", function(x) x@se_waic)

#' Create generic `se_waic<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric standard error of WAIC estimate to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' se_waic(m1) <- NaN
#' @export
setGeneric("se_waic<-", function(x, value) standardGeneric("se_waic<-"))

#' @rdname Model-class
setMethod("se_waic<-", "Model", function(x, value) {
  x@se_waic <- value
  x
})

#' Create generic `p_eff` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric estimate of the number of effective parameters when computing
#'   WAIC.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' p_eff(m1)
#' @export
setGeneric("p_eff", function(x) standardGeneric("p_eff"))

#' @rdname Model-class
setMethod("p_eff", "Model", function(x) x@p_eff)

#' Create generic `p_eff<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric value of effective number of parameters estimate from
#'   WAIC to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' p_eff(m1) <- NaN
#' @export
setGeneric("p_eff<-", function(x, value) standardGeneric("p_eff<-"))

#' @rdname Model-class
setMethod("p_eff<-", "Model", function(x, value) {
  x@p_eff <- value
  x
})

#' Create generic `lpd` function for class
#'
#' @param x An `Model` object.
#'
#' @return Numeric log-predictive density used in WAIC.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' lpd(m1)
#' @export
setGeneric("lpd", function(x) standardGeneric("lpd"))

#' @rdname Model-class
setMethod("lpd", "Model", function(x) x@lpd)

#' Create generic `lpd<-` function for class
#'
#' @param x An `Model` object.
#' @param value Numeric matrix of log predictive densities in each document to
#'   assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' lpd(m1) <- matrix(NaN, nrow = 1, ncol = 1)
#' @export
setGeneric("lpd<-", function(x, value) standardGeneric("lpd<-"))

#' @rdname Model-class
setMethod("lpd<-", "Model", function(x, value) {
  x@lpd <- value
  x
})

#' Create generic `extra` function for class
#'
#' @param x An `Model` object.
#'
#' @return A list of model fitting information including time elapsed, label
#'   switching correction status, and the original function call.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' extra(m1)
#' @export
setGeneric("extra", function(x) standardGeneric("extra"))

#' @rdname Model-class
setMethod("extra", "Model", function(x) x@extra)

#' Create generic `extra<-` function for class
#'
#' @param x An `Model` object.
#' @param value List of additional model fitting information to assign to slot.
#"
#' @return None.
setGeneric("extra<-", function(x, value) standardGeneric("extra<-"))

#' @rdname Model-class
setMethod("extra<-", "Model", function(x, value) {
  x@extra <- value
  x
})

#' Helper function (constructor) for `Model` class
#'
#' @return A [Model-class] object.
#'
#' @examples
#' m1 <- Model(ndocs = 1)
#' print(m1)
#'
#' @export
#' @rdname Model-class
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

# Validator function for `Model` class
setValidity("Model", function(object) {
  if ( (object@nchain != NROW(object@eta)) |
       (object@nchain != NROW(object@lpd))) {
    "@eta and @lpd must have the same number of rows (equal to @nchain)"
  } else if ( (object@nchain != NROW(object@eta)) |
              (object@nchain != length(object@loglike)) |
              (object@nchain != length(object@logpost)) |
              (object@nchain != NROW(object@lpd))) {
    "@loglike and @logpost must have the same number of rows (equal to @nchain)"
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

########----------------------------------------------------------------########
########                             Mlr Class                          ########
########----------------------------------------------------------------########

#' Create generic `sigma2` function for class
#'
#' @param x An `Mlr` object.
#'
#' @return Numeric vector of posterior draws of residual variance.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' sigma2(m1)
#' @export
setGeneric("sigma2", function(x) standardGeneric("sigma2"))

#' @rdname Mlr-class
setMethod("sigma2", "Mlr", function(x) x@sigma2)

#' Create generic `sigma2<-` function for class
#'
#' @param x An `Mlr` object.
#' @param value Numeric value of residual variance to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' sigma2(m1) <- 1.0
#' @export
setGeneric("sigma2<-", function(x, value) standardGeneric("sigma2<-"))

#' @rdname Mlr-class
setMethod("sigma2<-", "Mlr", function(x, value) {
  x@sigma2 <- value
  x
})

#' Create generic `a0` function for class
#'
#' @param x An `Mlr` object.
#'
#' @return Double value of shape prior parameter for residual variance.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' a0(m1)
#' @export
setGeneric("a0", function(x) standardGeneric("a0"))

#' @rdname Mlr-class
setMethod("a0", "Mlr", function(x) x@a0)

#' Create generic `a0<-` function for class
#'
#' @param x An `Mlr` object.
#' @param value Numeric shape parameter for residual variance prior to assign to
#'   slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' a0(m1) <- 1.0
#' @export
setGeneric("a0<-", function(x, value) standardGeneric("a0<-"))

#' @rdname Mlr-class
setMethod("a0<-", "Mlr", function(x, value) {
  x@a0 <- value
  x
})

#' Create generic `b0` function for class
#'
#' @param x An `Mlr` object.
#'
#' @return Double value of rate prior parameter for residual variance.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' b0(m1)
#' @export
setGeneric("b0", function(x) standardGeneric("b0"))

#' @rdname Mlr-class
setMethod("b0", "Mlr", function(x) x@b0)

#' Create generic `b0<-` function for class
#'
#' @param x An `Mlr` object.
#' @param value Numeric value of rate parameter for residual variance prior to
#'   assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' b0(m1) <- 1.0
#' @export
setGeneric("b0<-", function(x, value) standardGeneric("b0<-"))

#' @rdname Mlr-class
setMethod("b0<-", "Mlr", function(x, value) {
  x@b0 <- value
  x
})

#' Helper function (constructor) for `Mlr` class
#'
#' @return An [Mlr-class] object.
#'
#' @examples
#' m1 <- Mlr(ndocs = 1)
#' print(m1)
#'
#' @export
#' @rdname Mlr-class
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

########----------------------------------------------------------------########
########                          Logistic Class                        ########
########----------------------------------------------------------------########

#' Create generic `proposal_sd` function for class
#'
#' @param x An `Logistic` object.
#'
#' @return Numeric vector of proposal scales for Metropolis step for regression
#'   coefficients sampling.
#'
#' @examples
#' m1 <- Logistic(ndocs = 1)
#' proposal_sd(m1)
#' @export
setGeneric("proposal_sd", function(x) standardGeneric("proposal_sd"))

#' @rdname Logistic-class
setMethod("proposal_sd", "Logistic", function(x) x@proposal_sd)

#' Create generic `proposal_sd<-` function for class
#'
#' @param x An `Logistic` object.
#' @param value Numeric vector of scale parameters for Metropolis sampling of
#'   regression coefficients to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Logistic(ndocs = 1)
#' proposal_sd(m1) <- c(2.38, 2.38)
#' @export
setGeneric("proposal_sd<-", function(x, value) standardGeneric("proposal_sd<-"))

#' @rdname Logistic-class
setMethod("proposal_sd<-", "Logistic", function(x, value) {
  x@proposal_sd <- value
  validObject(x)
  x
})

#' Helper function (constructor) for `Logistic` class
#'
#' @return A [Logistic-class] object.
#'
#' @examples
#' m1 <- Logistic(ndocs = 1)
#' print(m1)
#'
#' @export
#' @rdname Logistic-class
Logistic <- function(proposal_sd = NaN, ...) {
  super <- Model(...)
  proposal_sd <- as.double(proposal_sd)

  new("Logistic", proposal_sd = proposal_sd,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

########----------------------------------------------------------------########
########                            Sldax Class                         ########
########----------------------------------------------------------------########

#' Create generic `topics` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Integer array of categorical topic labels for each word in each
#'   document across sampler iterations.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' topics(m1)
#' @export
setGeneric("topics", function(x) standardGeneric("topics"))

#' @rdname Sldax-class
setMethod("topics", "Sldax", function(x) x@topics)

#' Create generic `topics<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Integer array of topic assignment draws for each word to assign
#'   to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' topics(m1) <- array(c(2, 2, 2, 1), dim = c(1, 4, 1))
#' @export
setGeneric("topics<-", function(x, value) standardGeneric("topics<-"))

#' @rdname Sldax-class
setMethod("topics<-", "Sldax", function(x, value) {
  x@topics <- value
  x
})

#' Create generic `theta` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Numeric array of topic proportions for each document across sampler
#'   iterations.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' theta(m1)
#' @export
setGeneric("theta", function(x) standardGeneric("theta"))

#' @rdname Sldax-class
setMethod("theta", "Sldax", function(x) x@theta)

#' Create generic `theta<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Numeric array of topic proportions to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' theta(m1) <- array(c(0.5, 0.5), dim = c(1, 2, 1))
#' @export
setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))

#' @rdname Sldax-class
setMethod("theta<-", "Sldax", function(x, value) {
  x@theta <- value
  x
})

#' Create generic `beta_` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return A numeric array of topic-word probability distributions across
#'   sampler iterations.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' beta_(m1)
#' @export
setGeneric("beta_", function(x) standardGeneric("beta_"))

#' @rdname Sldax-class
setMethod("beta_", "Sldax", function(x) x@beta)

#' Create generic `beta_<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Numeric array of topic-word probabilities to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' beta_(m1) <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
#' @export
setGeneric("beta_<-", function(x, value) standardGeneric("beta_<-"))

#' @rdname Sldax-class
setMethod("beta_<-", "Sldax", function(x, value) {
  x@beta <- value
  x
})

#' Create generic `gamma_` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Double value of parameter for symmetric Dirichlet distribution prior
#'  on the topic-word probabilities.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' gamma_(m1)
#' @export
setGeneric("gamma_", function(x) standardGeneric("gamma_"))

#' @rdname Sldax-class
setMethod("gamma_", "Sldax", function(x) x@gamma)

#' Create generic `gamma_<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Numeric parameter for symmetric Dirichlet prior on topic-word
#'   probabilities to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' gamma_(m1) <- 1.0
#' @export
setGeneric("gamma_<-", function(x, value) standardGeneric("gamma_<-"))

#' @rdname Sldax-class
setMethod("gamma_<-", "Sldax", function(x, value) {
  x@gamma <- value
  x
})

#' Create generic `alpha` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Double value of parameter for symmetric Dirichlet distribution prior
#'  on the topic proportions.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' alpha(m1)
#' @export
setGeneric("alpha", function(x) standardGeneric("alpha"))

#' @rdname Sldax-class
setMethod("alpha", "Sldax", function(x) x@alpha)

#' Create generic `alpha<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Numeric parameter for symmetric Dirichlet prior on topic
#'   proportions to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' alpha(m1) <- 1.0
#' @export
setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))

#' @rdname Sldax-class
setMethod("alpha<-", "Sldax", function(x, value) {
  x@alpha <- value
  x
})

#' Create generic `ntopics` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Integer number of topics in model.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' ntopics(m1)
#' @export
setGeneric("ntopics", function(x) standardGeneric("ntopics"))

#' @rdname Sldax-class
setMethod("ntopics", "Sldax", function(x) x@ntopics)

#' Create generic `ntopics<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Integer number of topics to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' ntopics(m1) <- 2
#' @export
setGeneric("ntopics<-", function(x, value) standardGeneric("ntopics<-"))

#' @rdname Sldax-class
setMethod("ntopics<-", "Sldax", function(x, value) {
  x@ntopics <- value
  x
})

#' Create generic `nvocab` function for class
#'
#' @param x An `Sldax` object.
#'
#' @return Integer number of unique terms in vocabulary.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' nvocab(m1)
#' @export
setGeneric("nvocab", function(x) standardGeneric("nvocab"))

#' @rdname Sldax-class
setMethod("nvocab", "Sldax", function(x) x@nvocab)

#' Create generic `nvocab<-` function for class
#'
#' @param x An `Sldax` object.
#' @param value Numeric number of unique terms in vocabulary to assign to slot.
#'
#' @return None.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' nvocab(m1) <- 2L
#' @export
setGeneric("nvocab<-", function(x, value) standardGeneric("nvocab<-"))

#' @rdname Sldax-class
setMethod("nvocab<-", "Sldax", function(x, value) {
  x@nvocab <- value
  x
})

#' Helper function (constructor) for `Sldax` class
#'
#' @return A [Sldax-class] object.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' nvocab(m1) <- 2L
#' @export
#' @rdname Sldax-class
Sldax <- function(nvocab, topics, theta, beta, ntopics = 2.0, alpha = 1.0,
                  gamma = 1.0, ...) {
  super <- Model(...)
  nvocab <- as.double(nvocab)
  topics <- as.array(topics)
  theta <- as.array(theta)
  beta <- as.array(beta)
  ntopics <- as.double(ntopics)
  alpha <- as.double(alpha)
  gamma <- as.double(gamma)

  new("Sldax", nvocab = nvocab, topics = topics, theta = theta, beta = beta,
      ntopics = ntopics, alpha = alpha, gamma = gamma,
      ndocs = super@ndocs, nchain = super@nchain, mu0 = super@mu0,
      sigma0 = super@sigma0, eta_start = super@eta_start, eta = super@eta,
      loglike = super@loglike, logpost = super@logpost, waic = super@waic,
      se_waic = super@se_waic, p_eff = super@p_eff, lpd = super@lpd)
}

#' Summary functions for objects of class [Sldax-class]
#'
#' Obtain parameter estimates, model goodness-of-fit metrics, and posterior
#'   summaries.
#'
#' + `get_zbar()` computes empirical topic proportions from slot `@topics`.
#'
#' + `est_theta()` estimates the mean or median theta matrix.
#'
#' + `est_beta()` estimates the mean or median beta matrix.
#'
#' + `get_toptopics()` creates a [`tibble`][tibble::tibble()] of the topic
#'    proportion estimates for the top `ntopics` topics per document sorted by
#'    probability.
#'
#' + `get_topwords()` creates a [`tibble`][tibble::tibble()] of topics and the
#'    top `nwords` words per topic sorted by probability or term score.
#'
#' + `get_coherence()` computes the coherence metric for each topic (see Mimno,
#'    Wallach, Talley, Leenders, & McCallum, 2011).
#'
#' + `get_exclusivity()` computes the exclusivity metric for each topic (see
#'    Roberts, Stewart, & Airoldi, 2013).
#'
#' + `post_regression()` creates a [`coda::mcmc`][coda::mcmc] object
#' containing posterior information for the regression model parameters.
#'
#' + `gg_coef()` plots regression coefficients
#'
#'   + Warning: this function is deprecated.
#'
#'   + See `help("Deprecated")`.
#'
#' @name sldax-summary
NULL

#' Create generic `est_beta` function for class
#'
#' @param mcmc_fit An object of class [Sldax-class].
#' @param burn The number of draws to discard as a burn-in
#'   period (default: `0`).
#' @param thin The number of draws to skip as a thinning
#'   period (default: `1`; i.e., no thinning).
#' @param stat The summary statistic to use on the posterior
#'   draws (default: `"mean"`).
#'
#' @return A matrix of topic-word probability estimates.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
#' est_beta(m1, stat = "mean")
#' est_beta(m1, stat = "median")
#' @rdname sldax-summary
#' @export
setGeneric("est_beta",
           function(mcmc_fit,  burn = 0, thin = 1, stat = "mean") {

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if (length(dim(beta_(mcmc_fit))) != 3)
               stop("Only one draw of 'beta' available, so this function is not useful.")

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

             standardGeneric("est_beta")
           }
)

#' Create generic `est_theta` function for class
#'
#' @return A matrix of topic proportion estimates.
#'
#' @examples
#' m1 <- Sldax(ndocs = 2, nvocab = 2, nchain = 2,
#'             topics = array(c(1, 2, 2, 1,
#'                              1, 2, 2, 1), dim = c(2, 2, 2)),
#'             theta = array(c(0.5, 0.5,
#'                             0.5, 0.5,
#'                             0.5, 0.5,
#'                             0.5, 0.5), dim = c(2, 2, 2)),
#'             loglike = rep(NaN, times = 2),
#'             logpost = rep(NaN, times = 2),
#'             lpd = matrix(NaN, nrow = 2, ncol = 2),
#'             eta = matrix(0.0, nrow = 2, ncol = 2),
#'             mu0 = c(0.0, 0.0),
#'             sigma0 = diag(1, 2),
#'             eta_start = c(0.0, 0.0),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5,
#'                            0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2)))
#' est_theta(m1, stat = "mean")
#' est_theta(m1, stat = "median")
#' @rdname sldax-summary
#' @export
setGeneric("est_theta",
           function(mcmc_fit, burn = 0, thin = 1, stat = "mean") {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if (!is(mcmc_fit, "Sldax"))
               stop("'mcmc_fit' must be an Sldax object.")

             if (length(dim(theta(mcmc_fit))) != 3)
               stop("Only one draw of 'theta' available, so this function is not useful.")

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
             standardGeneric("est_theta")
           }
)

#' Create generic `get_coherence` function for class
#'
#' @param beta_ A \eqn{K} x \eqn{V} matrix of word-topic probabilities.
#'   Each row sums to 1.
#' @param docs The \eqn{D} x max(\eqn{N_d}) matrix of documents (word indices)
#'   used to fit the [Sldax-class] model.
#' @param nwords The number of highest-probability words per topic to consider
#'   where \eqn{M \le V} and \eqn{V} is the size of the corpus
#'   vocabulary. (default: `10`)
#'
#' @return A numeric vector of coherence scores for each
#'   topic (more positive is better).
#'
#' @examples
#' mdoc <- matrix(c(1, 2, 2, 1), nrow = 1)
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 2), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.4, 0.5, 0.6), dim = c(2, 2, 1)))
#' bhat <- est_beta(m1)
#' get_coherence(bhat, docs = mdoc, nwords = nvocab(m1))
#' @rdname sldax-summary
#' @export
setGeneric("get_coherence",
           function(beta_, docs, nwords = 10) {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(beta_))
               stop("Please supply an array to 'beta_'.")


             if ( length(dim(beta_)) != 2 )
               stop("'beta_' does not appear to be a K x V matrix.")

             if (any(beta_ < 0.0 | beta_ > 1.0))
               stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             tol <- 0.001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if ( !is.positive_integer(nwords) )
               stop("'nwords' must be a positive integer.")

             standardGeneric("get_coherence")
           }
)

#' Create generic `get_exclusivity` function for class
#'
#' @param weight The weight (between 0 and 1) to give to
#'   exclusivity (near 1) vs. frequency (near 0). (default: `0.7`).
#'
#' @return A numeric vector of exclusivity scores (more positive is better).
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 2), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.4, 0.5, 0.6), dim = c(2, 2, 1)))
#' bhat <- est_beta(m1)
#' get_exclusivity(bhat, nwords = nvocab(m1))
#' @rdname sldax-summary
#' @export
setGeneric("get_exclusivity",
           function(beta_, nwords = 10, weight = 0.7) {

             # passed_args <- names(as.list(match.call())[-1])

             if (missing(beta_))
               stop("Please supply an array to 'beta_'.")

             if (length(dim(beta_)) != 2)
               stop("'beta_' does not appear to be a K x V matrix.")

             if (any(beta_ < 0.0 | beta_ > 1.0))
               stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             tol <- 0.001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if ( !is.positive_integer(nwords) )
               stop("'nwords' must be a positive integer.")

             if ( ((weight >= 1.0) | (weight <= 0.0)))
               stop("'weight' must be between 0.0 and 1.0.")

             standardGeneric("get_exclusivity")
           }
)

#' Create generic `get_toptopics` function for class
#'
#' @param theta A D x K matrix of K topic proportions for all D documents.
#' @param ntopics The number of topics to retrieve (default: all topics).
#' @return A data frame of the `ntopics` most probable topics per document.
#'
#' @examples
#' m1 <- Sldax(ndocs = 2, nvocab = 2, nchain = 2,
#'             topics = array(c(1, 2, 2, 1,
#'                              1, 2, 2, 1), dim = c(2, 2, 2)),
#'             theta = array(c(0.4, 0.3,
#'                             0.6, 0.7,
#'                             0.45, 0.5,
#'                             0.55, 0.5), dim = c(2, 2, 2)),
#'             loglike = rep(NaN, times = 2),
#'             logpost = rep(NaN, times = 2),
#'             lpd = matrix(NaN, nrow = 2, ncol = 2),
#'             eta = matrix(0.0, nrow = 2, ncol = 2),
#'             mu0 = c(0.0, 0.0),
#'             sigma0 = diag(1, 2),
#'             eta_start = c(0.0, 0.0),
#'             beta = array(c(0.5, 0.5, 0.5, 0.5,
#'                            0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2)))
#' t_hat <- est_theta(m1, stat = "mean")
#' get_toptopics(t_hat, ntopics = ntopics(m1))
#' @rdname sldax-summary
#' @export
setGeneric("get_toptopics",
           function(theta, ntopics) {

             if (missing(theta)) stop("Please supply a matrix to 'theta'.")
             if (!is.matrix(theta)) stop("Please supply a matrix to 'theta'.")
             if (any(theta < 0.0 | theta > 1.0))
               stop("Entries of 'theta' must be between 0.0 and 1.0.")
             sum_rowsum_theta <- sum(rowSums(theta))
             d <- nrow(theta)
             K <- ncol(theta)
             tol <- 0.001
             if (sum_rowsum_theta > d + tol | sum_rowsum_theta < d - tol)
               stop("Rows of 'theta' must each sum to 1.0.")

             if (missing(ntopics)) ntopics <- K # Default
             if ( !is.positive_integer(ntopics) )
               stop("'ntopics' must be a positive integer.")

             standardGeneric("get_toptopics")
           }
)

#' Create generic `get_topwords` function for class
#'
#' @param nwords The number of words to retrieve (default: all).
#' @param vocab A character vector of length V containing the vocabulary.
#' @param method If `"termscore"`, use term scores (similar to tf-idf). If
#'  `"prob"`, use probabilities (default: `"termscore"`).
#'
#' @return A \eqn{K} x \eqn{V} matrix of term-scores (comparable to tf-idf).
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 2), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.4, 0.5, 0.6), dim = c(2, 2, 1)))
#' bhat <- est_beta(m1)
#' get_topwords(bhat, nwords = nvocab(m1), method = "termscore")
#' get_topwords(bhat, nwords = nvocab(m1), method = "prob")
#' @rdname sldax-summary
#' @export
setGeneric("get_topwords",
           function(beta_, nwords, vocab, method = "termscore") {

             if (missing(beta_)) stop("Argument 'beta_' is missing.")
             if (!is.matrix(beta_)) stop("Please supply a matrix to 'beta_'.")
             if (any(beta_ < 0.0 | beta_ > 1.0))
               stop("Entries of 'beta_' must be between 0.0 and 1.0.")
             sum_rowsum_beta <- sum(rowSums(beta_))
             K <- nrow(beta_)
             V <- ncol(beta_)
             tol <- .001
             if (sum_rowsum_beta > K + tol | sum_rowsum_beta < K - tol)
               stop("Rows of 'beta_' must each sum to 1.0.")

             if (missing(nwords)) nwords <- V # Default
             if ( !is.positive_integer(nwords) )
               stop("'nwords' must be a positive integer.")

             if (missing(vocab)) {
               vocab <- as.character(seq_len(V))
               message("'vocab' not supplied. Defaulting to
                        indices 1, 2, ..., V.")
             }
             if (length(vocab) < 2L)
               stop("'vocab' must contain at least two elements.")
             if (length(vocab) != V)
               stop("The number of elements in 'vocab' should equal the number of columns in 'beta_'.")
             if (!is.character(vocab))
               stop("'vocab' must be a character vector.")
             if (nwords > length(unique(vocab)))
               stop("'nwords' cannot exceed the number of unique terms in 'vocab'.")

             if (length(method) > 1) {
               method <- method[1]
               message("Multiple arguments were supplied to 'method'. Only using the first argument.")
             }
             if (!(method %in% c("termscore", "prob")))
               stop("'method' must be either 'termscore' or 'prob'.")

             standardGeneric("get_topwords")
           }
)

#' Create generic `get_zbar` function for class
#'
#' @return A matrix of empirical topic proportions per document.
#'
#' @examples
#' m1 <- Sldax(ndocs = 1, nvocab = 2,
#'             topics = array(c(1, 2, 2, 2), dim = c(1, 4, 1)),
#'             theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
#'             beta = array(c(0.5, 0.4, 0.5, 0.6), dim = c(2, 2, 1)))
#' get_zbar(m1)
#' @rdname sldax-summary
#' @export
setGeneric("get_zbar",
           function(mcmc_fit, burn = 0L, thin = 1L) {

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

             standardGeneric("get_zbar")
           }
)

#' Generic function to summarize regression relationships for `Mlr`,
#'   `Logistic`, or `Sldax` objects
#'
#' For SLDA or SLDAX models, label switching is handled during estimation in the
#' [`gibbs_sldax()`][gibbs_sldax()] function with argument `correct_ls`, so it
#' is not addressed by this function.
#'
#' @return An object of class [`coda::mcmc`][coda::mcmc] summarizing the
#'   posterior distribution of the regression coefficients and residual
#'   variance (if applicable). Convenience functions such
#'   as [`summary()`][coda::summary.mcmc()] and [`plot()`][coda::plot.mcmc()]
#'   can be used for posterior summarization.
#'
#' @examples
#' data(mtcars)
#' m1 <- gibbs_mlr(mpg ~ hp, data = mtcars, m = 2)
#' post_regression(m1)
#' @rdname sldax-summary
#' @export
setGeneric("post_regression",
           function(mcmc_fit) {

             if (missing(mcmc_fit))
               stop("Please supply an object to 'mcmc_fit'.")
             if ( !( is(mcmc_fit, "Sldax") |
                     is(mcmc_fit, "Mlr") |
                     is(mcmc_fit, "Logistic")) )
               stop("'mcmc_fit' must be an `Sldax` or `Mlr` or `Logistic` object.")

             if ( !is.null(extra(mcmc_fit)$call$model) ) {
               if (extra(mcmc_fit)$call$model == "lda")
                stop("The `lda` model does not contain regression parameters, so this function is not useful.")
             }

             if ( nrow(eta(mcmc_fit)) < 2L )
               stop("Only one draw of 'eta' available, so this function is not useful.")

             standardGeneric("post_regression")
           }
)
