context("beta parameter estimation")

test_that("est_beta() handles missing arguments", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles wrong model object class", {
  fit <- Mlr(ndocs = 1)
  expect_error(
    est_beta(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("est_beta() handles non-integer 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, burn = 1.1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles negative 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, burn = -1),
    "'burn' must be a non-negative integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'burn' longer than chain length", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("est_beta() handles non-integer 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(docs)),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, thin = 1.1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles negative 'thin'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, thin = -1),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit, thin = 0),
    "'thin' must be a positive integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'thin' longer than chain length minus 'burn'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1), dim = c(1, 4, 1))
  theta <- array(c(0.5, 0.5), dim = c(1, 2, 1))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1))
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               topics = topics, theta = theta, beta = beta_)
  expect_error(
    est_beta(mcmc_fit = fit, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles multiple values for 'stat'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(1, 4, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5), dim = c(1, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 1)
  loglike <- logpost <- rep(NaN, 2)
  mu0 <- c(0, 0)
  sigma0 <- diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               nchain = 2, topics = topics, theta = theta, beta = beta_,
               eta = eta, lpd = lpd, loglike = loglike, logpost = logpost,
               mu0 = mu0, sigma0 = sigma0, eta_start = eta_start)
  expect_message(
    est_beta(mcmc_fit = fit, stat = c("mean", "median")),
    "Multiple arguments were supplied to 'stat'. Only using the first argument.",
    fixed = TRUE
  )
})

test_that("est_beta() handles invalid 'stat'", {
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(1, 4, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5), dim = c(1, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 1)
  loglike <- logpost <- rep(NaN, 2)
  mu0 <- c(0, 0)
  sigma0 <- diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               nchain = 2, topics = topics, theta = theta, beta = beta_,
               eta = eta, lpd = lpd, loglike = loglike, logpost = logpost,
               mu0 = mu0, sigma0 = sigma0, eta_start = eta_start)
  expect_error(
    est_beta(mcmc_fit = fit, stat = 1),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit, stat = "maen"),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles missing 'stat' by defaulting to mean", {

  t1 <- t2 <- matrix(NaN, nrow = 2, ncol = 2)
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(1, 4, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5), dim = c(1, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 1)
  loglike <- logpost <- rep(NaN, 2)
  mu0 <- c(0, 0)
  sigma0 <- diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               nchain = 2, topics = topics, theta = theta, beta = beta_,
               eta = eta, lpd = lpd, loglike = loglike, logpost = logpost,
               mu0 = mu0, sigma0 = sigma0, eta_start = eta_start)
  for (iter in seq_len(nchain(fit))) {
    for (topic in seq_len(ntopics(fit))) {
      for (word in seq_len(nvocab(fit))) {
        if (iter == 1) t1[topic, word] <- sum(
          topics(fit)[, , iter] == topic & docs == word)
        if (iter == 2) t2[topic, word] <- sum(
          topics(fit)[, , iter] == topic & docs == word)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  expect_equal(
    est_beta(mcmc_fit = fit),
    apply(array(c(t1, t2), c(2, 2, 2)), c(1, 2), mean)
  )
})

test_that("est_beta() computes median estimate correctly", {
  t1 <- t2 <- matrix(NaN, nrow = 2, ncol = 2)
  docs <- matrix(c(1, 2, 1, 2), nrow = 1)
  topics <- array(c(1, 2, 2, 1,
                    1, 2, 2, 1), dim = c(1, 4, 2))
  theta <- array(c(0.5, 0.5,
                   0.5, 0.5), dim = c(1, 2, 2))
  beta_ <- array(c(0.5, 0.5, 0.5, 0.5,
                   0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 2))
  eta_start <- c(1, -1)
  eta <- matrix(c(1, -1, 1, -1), byrow = TRUE, nrow = 2)
  lpd <- matrix(NaN, nrow = 2, ncol = 1)
  loglike <- logpost <- rep(NaN, 2)
  mu0 <- c(0, 0)
  sigma0 <- diag(1, 2)
  fit <- Sldax(ndocs = nrow(docs), nvocab = length(unique(as.numeric(docs))),
               nchain = 2, topics = topics, theta = theta, beta = beta_,
               eta = eta, lpd = lpd, loglike = loglike, logpost = logpost,
               mu0 = mu0, sigma0 = sigma0, eta_start = eta_start)
  for (iter in seq_len(nchain(fit))) {
    for (topic in seq_len(ntopics(fit))) {
      for (word in seq_len(nvocab(fit))) {
        if (iter == 1) t1[topic, word] <- sum(
          topics(fit)[, , iter] == topic & docs == word)
        if (iter == 2) t2[topic, word] <- sum(
          topics(fit)[, , iter] == topic & docs == word)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  expect_equal(
    est_beta(mcmc_fit = fit, stat = "median"),
    apply(array(c(t1, t2), c(2, 2, 2)), c(1, 2), median)
  )
})
