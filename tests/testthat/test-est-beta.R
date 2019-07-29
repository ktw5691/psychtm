context("beta parameter estimation")

test_that("est_beta() handles missing arguments", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2; fit@nvocab <- 2
  expect_error(
    est_beta(docs = docs),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit),
    "Please supply a matrix for 'docs'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles wrong model object class", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Mlr()
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'docs' with wrong number of rows", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 3L
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs),
    "'docs' must have 3 rows.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'docs' with wrong number of columns", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@topics <- array(c(1, 2, 1, 2, 1, 0, 0, 0,
                        1, 2, 1, 2, 1, 0, 0, 0),
                      dim = c(2, 4, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs),
    "'docs' must have 4 columns.",
    fixed = TRUE
  )
})

test_that("est_beta() handles non-integer 'burn'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, burn = 1.1),
    "'burn' must be an integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles negative 'burn'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, burn = -1),
    "'burn' must be non-negative.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'burn' longer than chain length", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@nchain <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("est_beta() handles non-integer 'thin'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, thin = 1.1),
    "'thin' must be an integer.",
    fixed = TRUE
  )
})

test_that("est_beta() handles negative 'thin'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, thin = -1),
    "'thin' must be positive.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, thin = 0),
    "'thin' must be positive.",
    fixed = TRUE
  )
})

test_that("est_beta() handles 'thin' longer than chain length minus 'burn'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs <- 2L
  fit@nchain <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, burn = 1, thin = 1),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, burn = 1, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles multiple values for 'stat'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 2L
  fit@nvocab <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_message(
    est_beta(mcmc_fit = fit, docs = docs, stat = c("mean", "median")),
    "Multiple arguments were supplied to 'stat'. Only using the first argument.",
    fixed = TRUE
  )
})

test_that("est_beta() handles invalid 'stat'", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 2L
  fit@nvocab <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, stat = 1),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
  expect_error(
    est_beta(mcmc_fit = fit, docs = docs, stat = "maen"),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
})

test_that("est_beta() handles missing 'stat' by defaulting to mean", {
  docs <- matrix(c(1, 2, 0,
                   2, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@nvocab <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0,
                        2, 1, 2, 0, 0, 0),
                      dim = c(2, 3, 3))
  fit@gamma <- 0.0
  t1 <- t2 <- t3 <- matrix(NaN, nrow = 2, ncol = 2)
  for (iter in seq_len(fit@nchain)) {
    for (topic in seq_len(fit@ntopics)) {
      for (word in seq_len(fit@nvocab)) {
        if (iter == 1) t1[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
        if (iter == 2) t2[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
        if (iter == 3) t3[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  t3 <- t3 / rowSums(t3)
  expect_equal(
    est_beta(mcmc_fit = fit, docs = docs),
    apply(array(c(t1, t2, t3), c(2, 2, 3)), c(1, 2), mean)
  )
})

test_that("est_beta() computes median estimate correctly", {
  docs <- matrix(c(1, 2, 0,
                   1, 0, 0),
                 nrow = 2, byrow = TRUE)
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@nvocab <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0,
                        2, 1, 2, 0, 0, 0),
                      dim = c(2, 3, 3))
  fit@gamma <- 0.0
  t1 <- t2 <- t3 <- matrix(NaN, nrow = 2, ncol = 2)
  for (iter in seq_len(fit@nchain)) {
    for (topic in seq_len(fit@ntopics)) {
      for (word in seq_len(fit@nvocab)) {
        if (iter == 1) t1[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
        if (iter == 2) t2[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
        if (iter == 3) t3[topic, word] <- sum(fit@topics[, , iter] == topic & docs == word)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  t3 <- t3 / rowSums(t3)
  expect_equal(
    est_beta(mcmc_fit = fit, docs = docs, stat = "median"),
    apply(array(c(t1, t2, t3), c(2, 2, 3)), c(1, 2), median)
  )
})
