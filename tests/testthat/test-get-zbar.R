context("empirical topic proportion estimation")

test_that("get_zbar() handles missing 'mcmc_fit'", {
  expect_error(
    get_zbar(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles wrong model object class", {
  fit <- Mlr(ndocs = 1)
  expect_error(
    get_zbar(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles non-integer 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 1.1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles negative 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, burn = -1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles 'burn' longer than chain length", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles non-integer 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, thin = 1.1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles negative 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, thin = -1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
  expect_error(
    get_zbar(mcmc_fit = fit, thin = 0),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles 'thin' longer than chain length minus 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    get_zbar(mcmc_fit = fit, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("get_zbar() correctly computes empirical topic proportions", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 1, 1, 2), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)

  zbar <- c(0.75, 0.25)
  expect_equal(
    get_zbar(mcmc_fit = fit),
    zbar
  )
})
