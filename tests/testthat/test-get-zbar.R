context("empirical topic proportion estimation")

test_that("get_zbar() handles missing 'mcmc_fit'", {
  expect_error(
    get_zbar(),
    "Please supply an object to 'mcmc_fit'.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles wrong model object class", {
  fit <- Mlr()
  expect_error(
    get_zbar(mcmc_fit = fit),
    "'mcmc_fit' must be an Sldax object.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles non-integer 'burn'", {
  fit <- Sldax()
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 1.1),
    "'burn' must be an integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles negative 'burn'", {
  fit <- Sldax()
  expect_error(
    get_zbar(mcmc_fit = fit, burn = -1),
    "'burn' must be non-negative.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles 'burn' longer than chain length", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 3),
    "'burn' cannot exceed length of chain.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles non-integer 'thin'", {
  fit <- Sldax()
  expect_error(
    get_zbar(mcmc_fit = fit, thin = 1.1),
    "'thin' must be an integer.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles negative 'thin'", {
  fit <- Sldax()
  expect_error(
    get_zbar(mcmc_fit = fit, thin = -1),
    "'thin' must be positive.",
    fixed = TRUE
  )
  expect_error(
    get_zbar(mcmc_fit = fit, thin = 0),
    "'thin' must be positive.",
    fixed = TRUE
  )
})

test_that("get_zbar() handles 'thin' longer than chain length minus 'burn'", {
  fit <- Sldax()
  fit@nchain <- 2L
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 1, thin = 1),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
  expect_error(
    get_zbar(mcmc_fit = fit, burn = 1, thin = 2),
    "'thin' cannot exceed length of chain less 'burn'.",
    fixed = TRUE
  )
})

test_that("get_zbar() correctly computes empirical topic proportions", {
  fit <- Sldax()
  fit@ndocs  <- 2L
  fit@nchain <- 3L
  fit@topics <- array(c(1, 2, 1, 0, 0, 0,
                        1, 2, 1, 0, 0, 0,
                        2, 1, 2, 0, 0, 0),
                      dim = c(2, 3, 3))
  zbar <- matrix(c(1.0, 0.0,
                   0.0, 1.0), nrow = 2, byrow = TRUE)
  expect_equal(
    get_zbar(mcmc_fit = fit),
    zbar
  )
})
