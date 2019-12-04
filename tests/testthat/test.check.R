# @author Quentin Grimonprez

context("Check functions")

test_that("checkData works", {
  
  expect_error(checkData(c()), "data must be a data.frame.")
  
  df <- data.frame(id = 1)
  expect_error(checkData(df), "Missing columns in data: time, state.")
  
  df <- data.frame(id = 1, time = 1, state = "e")
  expect_error(checkData(df), "There is only one row or less.")
  
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
  expect_equal(is.whole.number(NA), NA)
  expect_true(all(is.whole.number(1:5)))
  expect_equal(is.whole.number(c(0.5, 2)), c(FALSE, TRUE))
})
  

test_that("checkDataEndTmax works", {
  df <- data.frame(id = 1:2, time = 1:2, state = 1:2)
  expect_error(checkDataEndTmax(df), "Each individual must end with the same time value.")
  
  df <- data.frame(id = 1:2, time = rep(1, 2), state = 1:2)
  expect_silent(checkDataEndTmax(df))
})

test_that("checkLogical works", {
  expect_silent(checkLogical(TRUE, "aa"))
  expect_error(checkLogical(c(TRUE, TRUE), "aa"), "aa must be either TRUE or FALSE.")
  expect_error(checkLogical(3, "aa"), "aa must be either TRUE or FALSE.")
})
