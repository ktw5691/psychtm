context("plotting of regression coefficients")

test_that("gg_coef() handles missing 'mcmc_fit'", {
  expect_error(
    gg_coef(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles wrong model object class", {
  fit <- Mlr(ndocs = 1)
  expect_error(
    gg_coef(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles non-integer 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 1.1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles negative 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, burn = -1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles 'burn' longer than chain length", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles non-integer 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, thin = 1.1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles negative 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, thin = -1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
  expect_error(
    gg_coef(mcmc_fit = fit, thin = 0),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles 'thin' longer than chain length minus 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles multiple values for 'stat'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 2)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(2, 2, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5), dim = c(2, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 2)
  loglike <- logpost <- rep(NaN, 2)
  mu0 = c(0, 0)
  sigma0 = diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))), nchain = 2,
               topics = topics, theta = theta, beta = beta_, eta = eta,
               lpd = lpd, loglike = loglike, logpost = logpost, mu0 = mu0,
               sigma0 = sigma0, eta_start = eta_start)
  colnames(eta(fit)) <- c("topic1", "topic2")
  expect_message(
    gg_coef(mcmc_fit = fit, stat = c("mean", "median")),
    "Multiple arguments were supplied to 'stat'. Only using the first argument.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles invalid 'stat'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(mcmc_fit = fit, stat = 1),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
  expect_error(
    gg_coef(mcmc_fit = fit, stat = "maen"),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles missing 'stat' by defaulting to mean", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 2)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(2, 2, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5), dim = c(2, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 2)
  loglike <- logpost <- rep(NaN, 2)
  mu0 = c(0, 0)
  sigma0 = diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))), nchain = 2,
               topics = topics, theta = theta, beta = beta_, eta = eta,
               lpd = lpd, loglike = loglike, logpost = logpost, mu0 = mu0,
               sigma0 = sigma0, eta_start = eta_start)
  colnames(eta(fit)) <- c("topic1", "topic2")
  expect_equal(
    gg_coef(mcmc_fit = fit)$data$est,
    c(topic2 = -1, topic1 = 1) # Sorted by coefficient rank
  )
})

test_that("gg_coef() correctly computes median estimate of eta", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 2)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(2, 2, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5,
                   0.5, 0.5), dim = c(2, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 2)
  loglike <- logpost <- rep(NaN, 2)
  mu0 = c(0, 0)
  sigma0 = diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))), nchain = 2,
               topics = topics, theta = theta, beta = beta_, eta = eta,
               lpd = lpd, loglike = loglike, logpost = logpost, mu0 = mu0,
               sigma0 = sigma0, eta_start = eta_start)
  colnames(eta(fit)) <- c("topic1", "topic2")
  expect_equal(
    gg_coef(mcmc_fit = fit, stat = "median")$data$est,
    c(topic2 = -1, topic1 = 1) # Sorted by coefficient rank
  )
})

test_that("'errorbw' is numeric", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(fit, errorbw = "a"),
    "'errorbw' must be numeric."
  )
})

test_that("'errorbw' is positive", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    gg_coef(fit, errorbw = -1),
    "'errorbw' must be positive."
  )
})
