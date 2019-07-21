context("argument checking helper functions")

test_that("check_int() handles integers", {
  expect_true(
    check_int(1L) || check_int(1.0) || check_int(-1L) || check_int(-1.0)
  )
})

test_that("check_int() handles non-integers", {
  expect_false(
    check_int("a") ||
      check_int(1.1) ||
      check_int(-1.1) ||
      check_int(Inf) ||
      check_int(NaN) ||
      check_int(NA) ||
      check_int(TRUE)
  )
})

test_that("check_logical() handles booleans", {
  expect_true(
    check_logical(TRUE) || check_logical(FALSE)
  )
})

test_that("check_logical() handles non-booleans", {
  expect_false(
    check_logical("a") ||
      check_logical(1L) ||
      check_logical(1.0) ||
      check_logical(1.1) ||
      check_logical(-1.1) ||
      check_logical(Inf) ||
      check_logical(NaN) ||
      check_logical(NA)
  )
})
