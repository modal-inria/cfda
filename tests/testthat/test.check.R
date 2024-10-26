# @author Quentin Grimonprez

context("Check functions")

test_that("checkData works", {
  expect_error(checkData(c()), "data must be a data.frame.")

  df <- data.frame(id = 1)
  expect_error(checkData(df), "Missing columns in data: time, state.")

  df <- data.frame(id = 1, time = 1, state = "e")
  expect_error(checkData(df), "There is 1 row or less.")

  df <- data.frame(id = 1, time = 1, state = "e")
  expect_error(checkData(df, minSize = 10), "There is 10 row or less.")

  df <- data.frame(id = 1:2, time = c(1, NA), state = 1:2)
  expect_error(checkData(df), "There is some missing values.")

  df <- data.frame(id = 1:2, time = c(1, NA), state = letters[1:2])
  expect_error(checkData(df), "There is some missing values.")

  df <- data.frame(id = 1:2, time = 1:2, state = 1:2)
  expect_silent(checkData(df))
})


test_that("is.whole.number works", {
  expect_true(is.whole.number(2))
  expect_true(is.whole.number(2.))
  expect_false(is.whole.number(2.5))
  expect_equal(is.whole.number(NA), FALSE)
  expect_equal(is.whole.number("aaa"), FALSE)
  expect_true(all(is.whole.number(1:5)))
  expect_equal(is.whole.number(c(0.5, 2)), c(FALSE, TRUE))
})


test_that("checkDataEndTmax works", {
  df <- data.frame(id = 1:2, time = 1:2, state = 1:2)
  expect_error(checkDataEndTmax(df), "Each individual must end with the same time value.")

  df <- data.frame(id = 1:2, time = rep(1, 2), state = 1:2)
  expect_silent(checkDataEndTmax(df))
})

test_that("checkDataBeginTime works", {
  df <- data.frame(id = 1:2, time = 1:2, state = 1:2)
  expect_error(checkDataBeginTime(df), "Each individual must begin with the same time value.")

  df <- data.frame(id = 1:2, time = rep(1, 1), state = 1:2)
  expect_silent(checkDataBeginTime(df))
})

test_that("checkDataNoDuplicatedTimes works", {
  df <- data.frame(id = rep(1:2, each = 2), time = c(0, 0, 1, 2), state = rep(1:2, 2))
  expect_warning(checkDataNoDuplicatedTimes(df), "Some ids contain duplicated time values.")

  df <- data.frame(id = rep(1:2, each = 2), time = rep(1:2, 2), state = rep(1:2, 2))
  expect_silent(checkDataNoDuplicatedTimes(df))
})

test_that("checkLogical works", {
  expect_silent(checkLogical(TRUE, "aa"))
  expect_error(checkLogical(c(TRUE, TRUE), "aa"), "aa must be either TRUE or FALSE.")
  expect_error(checkLogical(3, "aa"), "aa must be either TRUE or FALSE.")
  expect_error(checkLogical(NA, "aa"), "aa must be either TRUE or FALSE.")
  expect_error(checkLogical(NaN, "aa"), "aa must be either TRUE or FALSE.")
})

test_that("checkInteger accepts valid integers within bounds", {
  expect_silent(checkInteger(5, minValue = 1, maxValue = 10))
  expect_silent(checkInteger(-3, minValue = -5, maxValue = 0))
  expect_silent(checkInteger(0))
})

test_that("checkInteger throws an error for non-integer values", {
  expect_error(checkInteger(5.5), "x must be an integer")
  expect_error(checkInteger("text", paramName = "a"), "a must be an integer")
  expect_error(checkInteger("text", minValue = 0, paramName = "a"), "a must be an integer")
})

test_that("checkInteger accepts NULL when acceptNULL is TRUE", {
  expect_silent(checkInteger(NULL, acceptNULL = TRUE))
})

test_that("checkInteger throws an error for values out of bounds", {
  expect_error(
    checkInteger(11, minValue = 1, maxValue = 10, minEqual = FALSE, maxEqual = FALSE),
    "x must be an integer > 1 and < 10"
  )
  expect_error(
    checkInteger(10, minValue = 1, maxValue = 10, minEqual = FALSE, maxEqual = FALSE),
    "x must be an integer > 1 and < 10"
  )
  expect_error(
    checkInteger(0, minValue = 1, maxValue = 10, minEqual = FALSE, maxEqual = FALSE),
    "x must be an integer > 1 and < 10"
  )
  expect_error(
    checkInteger(1, minValue = 1, maxValue = 10, minEqual = FALSE, maxEqual = FALSE),
    "x must be an integer > 1 and < 10"
  )
  expect_error(
    checkInteger(-6, minValue = -5, maxValue = 0, minEqual = FALSE, maxEqual = FALSE),
    "x must be an integer > -5 and < 0"
  )
})

test_that("checkInteger throws an error for values out of bounds (minEqual = TRUE, maxEqual = TRUE)", {
  expect_error(
    checkInteger(11, minValue = 1, maxValue = 10, minEqual = TRUE, maxEqual = TRUE),
    "x must be an integer >= 1 and <= 10"
  )
  expect_silent(checkInteger(10, minValue = 1, maxValue = 10, minEqual = FALSE, maxEqual = TRUE))
  expect_error(
    checkInteger(0, minValue = 1, maxValue = 10, minEqual = TRUE, maxEqual = TRUE),
    "x must be an integer >= 1 and <= 10"
  )
  expect_silent(checkInteger(1, minValue = 1, maxValue = 10, minEqual = TRUE, maxEqual = FALSE))
  expect_error(
    checkInteger(-6, minValue = -5, maxValue = 0, minEqual = TRUE, maxEqual = TRUE),
    "x must be an integer >= -5 and <= 0"
  )
})

test_that("checkInteger throws an error when wrong type", {
  expect_error(checkInteger(NULL), "x must be an integer")
  expect_error(checkInteger(NA), "x must be an integer")
  expect_error(checkInteger(c(1, 2, 3)), "x must be an integer")
})

test_that("checkInteger throws custom message", {
  expect_error(checkInteger(NULL, customMessage = "aa"), "aa")
})
