context("argument checking helper functions")

test_that("is.whole_number() handles integers", {
  expect_true(
    is.whole_number(1L) &
      is.whole_number(1.0) & is.whole_number(-1L) &
      is.whole_number(-1.0)
  )
})

test_that("is.whole_number() handles non-integers", {
  expect_false(
      is.whole_number(1.1) |
      is.whole_number(-1.1)
  )
})

test_that("is.positive_integer() handles integers", {
  expect_true(
    is.positive_integer(1L) &
      is.positive_integer(1.0)
  )
})

test_that("is.positive_integer() handles non-positive integers", {
  expect_false(
    is.positive_integer(0.0) |
      is.positive_integer(0L) |
      is.positive_integer(-1.0) |
      is.positive_integer(-1L)
  )
})

test_that("is.non_negative_integer() handles non-negative integers", {
  expect_true(
    is.non_negative_integer(0L) &
      is.non_negative_integer(0.0) &
      is.non_negative_integer(1L) &
      is.non_negative_integer(1.0)
  )
})

test_that("is.non_negative_integer() handles negative integers", {
  expect_false(
      is.non_negative_integer(-1.0) |
      is.non_negative_integer(-1L)
  )
})

test_that("check_logical() handles booleans", {
  expect_true(
    check_logical(TRUE) & check_logical(FALSE)
  )
})

test_that("check_logical() handles non-booleans", {
  expect_false(
    check_logical("a") |
      check_logical(1L) |
      check_logical(1.0) |
      check_logical(1.1) |
      check_logical(-1.1) |
      check_logical(Inf) |
      check_logical(NaN) |
      check_logical(NA)
  )
})
