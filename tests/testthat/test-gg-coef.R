context("plotting of regression coefficients")

test_that("gg_coef() handles missing 'mcmc_fit'", {
  expect_error(
    gg_coef(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles wrong model object class", {
  fit <- Mlr()
  expect_error(
    gg_coef(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles non-integer 'burn'", {
  fit <- Sldax()
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 1.1),
    "'burn' must be an integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles negative 'burn'", {
  fit <- Sldax()
  expect_error(
    gg_coef(mcmc_fit = fit, burn = -1),
    "'burn' must be non-negative.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles 'burn' longer than chain length", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles non-integer 'thin'", {
  fit <- Sldax()
  expect_error(
    gg_coef(mcmc_fit = fit, thin = 1.1),
    "'thin' must be an integer.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles negative 'thin'", {
  fit <- Sldax()
  expect_error(
    gg_coef(mcmc_fit = fit, thin = -1),
    "'thin' must be positive.",
    fixed = TRUE
  )
  expect_error(
    gg_coef(mcmc_fit = fit, thin = 0),
    "'thin' must be positive.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles 'thin' longer than chain length minus 'burn'", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 1, thin = 1),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
  expect_error(
    gg_coef(mcmc_fit = fit, burn = 1, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles multiple values for 'stat'", {
  fit <- Sldax()
  fit@nchain <- 2L
  fit@ndocs  <- 2L
  fit@eta <- matrix(c(0, 0,
                      0.1, 0.1),
                    nrow = 2, byrow = TRUE)
  colnames(fit@eta) <- c("topic1", "topic2")
  expect_message(
    gg_coef(mcmc_fit = fit, stat = c("mean", "median")),
    "Multiple arguments were supplied to 'stat'. Only using the first argument.",
    fixed = TRUE
  )
})

test_that("gg_coef() handles invalid 'stat'", {
  fit <- Sldax()
  fit@nchain <- 2L
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
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@eta <- matrix(c(0, 0,
                      0, 0.1,
                      0.1, 0),
                    nrow = 3, byrow = TRUE)
  colnames(fit@eta) <- c("topic1", "topic2")
  expect_equal(
    gg_coef(mcmc_fit = fit)$data$est,
    c(0.1 / 3, 0.1 / 3)
  )
})

test_that("gg_coef() correctly computes median estimate of eta", {
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@eta <- matrix(c(0, 0,
                      0, 0.1,
                      0.1, 0),
                    nrow = 3, byrow = TRUE)
  colnames(fit@eta) <- c("topic1", "topic2")
  expect_equal(
    gg_coef(mcmc_fit = fit, stat = "median")$data$est,
    c(0, 0)
  )
})
