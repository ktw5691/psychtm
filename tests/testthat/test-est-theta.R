context("theta parameter estimation")

test_that("est_theta() handles missing 'mcmc_fit'", {
  expect_error(
    est_theta(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("est_theta() handles wrong model object class", {
  fit <- Mlr()
  expect_error(
    est_theta(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("est_theta() handles non-integer 'burn'", {
  fit <- Sldax()
  expect_error(
    est_theta(mcmc_fit = fit, burn = 1.1),
    "'burn' must be an integer.",
    fixed = TRUE
  )
})

test_that("est_theta() handles negative 'burn'", {
  fit <- Sldax()
  expect_error(
    est_theta(mcmc_fit = fit, burn = -1),
    "'burn' must be non-negative.",
    fixed = TRUE
  )
})

test_that("est_theta() handles 'burn' longer than chain length", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    est_theta(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("est_theta() handles non-integer 'thin'", {
  fit <- Sldax()
  expect_error(
    est_theta(mcmc_fit = fit, thin = 1.1),
    "'thin' must be an integer.",
    fixed = TRUE
  )
})

test_that("est_theta() handles negative 'thin'", {
  fit <- Sldax()
  expect_error(
    est_theta(mcmc_fit = fit, thin = -1),
    "'thin' must be positive.",
    fixed = TRUE
  )
  expect_error(
    est_theta(mcmc_fit = fit, thin = 0),
    "'thin' must be positive.",
    fixed = TRUE
  )
})

test_that("est_theta() handles 'thin' longer than chain length minus 'burn'", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    est_theta(mcmc_fit = fit, burn = 1, thin = 1),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
  expect_error(
    est_theta(mcmc_fit = fit, burn = 1, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("est_theta() handles multiple values for 'stat'", {
  fit <- Sldax()
  fit@nchain <- 2L
  fit@ndocs  <- 2L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0),
                      dim = c(2, 3, 2))
  expect_message(
    est_theta(mcmc_fit = fit, stat = c("mean", "median")),
    "Multiple arguments were supplied to 'stat'. Only using the first argument.",
    fixed = TRUE
  )
})

test_that("est_theta() handles invalid 'stat'", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    est_theta(mcmc_fit = fit, stat = 1),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
  expect_error(
    est_theta(mcmc_fit = fit, stat = "maen"),
    "'stat' must be either 'mean' or 'median'.",
    fixed = TRUE
  )
})

test_that("est_theta() handles missing 'stat' by defaulting to mean", {
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0,
                        2, 1, 2, 0, 0, 0),
                      dim = c(2, 3, 3))
  fit@alpha <- 0.0
  t1 <- t2 <- t3 <- matrix(NaN, nrow = 2, ncol = 2)
  for (iter in seq_len(fit@nchain)) {
    for (topic in seq_len(fit@ntopics)) {
      for (doc in seq_len(fit@ndocs)) {
        if (iter == 1) t1[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
        if (iter == 2) t2[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
        if (iter == 3) t3[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  t3 <- t3 / rowSums(t3)
  expect_equal(
    est_theta(mcmc_fit = fit),
    apply(array(c(t1, t2, t3), c(2, 2, 3)), c(1, 2), mean)
  )
})

test_that("est_theta() correctly computes median", {
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0,
                        2, 1, 2, 0, 0, 0),
                      dim = c(2, 3, 3))
  fit@alpha <- 0.0
  t1 <- t2 <- t3 <- matrix(NaN, nrow = 2, ncol = 2)
  for (iter in seq_len(fit@nchain)) {
    for (topic in seq_len(fit@ntopics)) {
      for (doc in seq_len(fit@ndocs)) {
        if (iter == 1) t1[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
        if (iter == 2) t2[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
        if (iter == 3) t3[doc, topic] <- sum(fit@topics[doc, , iter] == topic)
      }
    }
  }
  t1 <- t1 / rowSums(t1)
  t2 <- t2 / rowSums(t2)
  t3 <- t3 / rowSums(t3)
  expect_equal(
    est_theta(mcmc_fit = fit, stat = "median"),
    apply(array(c(t1, t2, t3), c(2, 2, 3)), c(1, 2), median)
  )
})
