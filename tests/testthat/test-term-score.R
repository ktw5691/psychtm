context("term score calculation from beta")

test_that("term_score() handles missing argument", {
  expect_error(
    term_score(),
    "Please supply a matrix to 'beta_'.",
    fixed = TRUE
  )
})

test_that("term_score() handles wrong class", {
  fit <- data.frame()
  expect_error(
    term_score(fit),
    "Please supply a matrix to 'beta_'.",
    fixed = TRUE
  )
})

test_that("term_score() handles non-probability vector rows", {
  fit <- matrix(c(1, 0.1,
                  0.5, 0.5),
                nrow = 2, byrow = TRUE)
  expect_error(
    term_score(fit),
    "Rows of 'beta_' must each sum to 1.0.",
    fixed = TRUE
  )
  fit <- matrix(c(1.1, 0,
                  0.5, 0.5),
                nrow = 2, byrow = TRUE)
  expect_error(
    term_score(fit),
    "Entries of 'beta_' must be between 0.0 and 1.0.",
    fixed = TRUE
  )
})

test_that("term_score() is computed correctly", {
  mbeta <- matrix(c(0.8, 0.2,
                    0.5, 0.5),
                  nrow = 2, byrow = TRUE)
  denom <- apply(mbeta, 2, prod) ^ 0.5 # Product over topics (rows)
  mdenom <- rbind(denom, denom, deparse.level = 0)
  ts <- mbeta * log(mbeta / mdenom)
  expect_equal(
    term_score(mbeta),
    ts
  )
  mbeta <- matrix(c(0.5, 0.5,
                    0.5, 0.5),
                  nrow = 2)
  expect_equal(
    term_score(mbeta),
    matrix(0, nrow = 2, ncol = 2)
  )
})
