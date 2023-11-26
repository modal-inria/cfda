# Author: Quentin Grimonprez

context("Format Multivariate Data")

test_that("convert2mvcfd works with 1 indiv and different time", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 0.9), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(1, 5), time = c(0, 0.25, 0.5, 0.9, 1), state1 = c(1, 1, 2, 2, 1), state2 = c(1, 2, 2, 3, 3)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with several ind", {
  x1 <- data.frame(id = c(1, 1, 1, 2, 2), time = c(0, 0.5, 1, 0, 1.5), state = c(1, 2, 1, 1, 2))
  x2 <- data.frame(id = c(1, 1, 1, 2, 2), time = c(0, 0.25, 0.9, 0, 2), state = c(1, 2, 3, 1, 2))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(c(1, 2), c(5, 3)),
    time = c(0, 0.25, 0.5, 0.9, 1, 0, 1.5, 2),
    state1 = c(1, 1, 2, 2, 1, 1, 2, 2),
    state2 = c(1, 2, 2, 3, 3, 1, 1, 2)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with 1 indiv and same time", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 1), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(id = c(1, 1, 1, 1), time = c(0, 0.25, 0.5, 1), state1 = c(1, 1, 2, 1), state2 = c(1, 2, 2, 3))
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with more than 2 dataframes", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 1), state = c(1, 2, 3))
  x3 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 3, 1))
  x <- list(x1, x2, x3)

  expectedOut <- data.frame(
    id = c(1, 1, 1, 1), time = c(0, 0.25, 0.5, 1), state1 = c(1, 1, 2, 1), state2 = c(1, 2, 2, 3), state3 = c(1, 1, 3, 1)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd throws error with bad input", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))

  expect_error(convert2mvcfd(x1), regexp = "data must be a list of data.frames")
})

test_that("convert2mvcfd works with stateColumns", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 1), state = c(1, 2, 3))
  x3 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 3, 1))
  x <- list(x1, x2, x3)

  expectedOut <- data.frame(
    id = c(1, 1, 1, 1), time = c(0, 0.25, 0.5, 1), dimA = c(1, 1, 2, 1), dimB = c(1, 2, 2, 3), dimC = c(1, 1, 3, 1)
  )
  out <- convert2mvcfd(x, stateColumns = c("dimA", "dimB", "dimC"))
  expect_equal(out, expectedOut)
})

test_that("convertMvcfd2cfd works", {
  # 2 columns
  x <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state1 = c(1, 2, 1), state2 = c(2, 1, 3))
  expectedOut <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = as.factor(c("1_2", "2_1", "1_3")))

  out <- convertMvcfd2cfd(x)
  expect_equal(out, expectedOut)

  # > 2 columns
  x <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state1 = c(1, 2, 1), state2 = c(2, 1, 3), state3 = c("a", "b", "c"))
  expectedOut <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = as.factor(c("1_2_a", "2_1_b", "1_3_c")))

  out <- convertMvcfd2cfd(x)
  expect_equal(out, expectedOut)
})

test_that("convertMvcfd2cfd works with custom sep", {
  x <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state1 = c(1, 2, 1), state2 = c(2, 1, 3))
  expectedOut <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = as.factor(c("1 & 2", "2 & 1", "1 & 3")))

  out <- convertMvcfd2cfd(x, sep = " & ")
  expect_equal(out, expectedOut)
})

test_that("convertMvcfd2cfd works with custom stateColumns", {
  x <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state1 = c(1, 2, 1), state2 = c(2, 1, 3))
  expectedOut <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = as.factor(c("1 & 2", "2 & 1", "1 & 3")))

  out <- convertMvcfd2cfd(x, sep = " & ")
  expect_equal(out, expectedOut)
})

test_that("convertMvcfd2cfd throws error with bad input", {
  x1 <- data.frame(id1 = c(1, 1, 1), time = c(0, 0.5, 1))

  expect_error(convertMvcfd2cfd(2), regexp = "data must be a data.frame")
  expect_error(convertMvcfd2cfd(x1), regexp = "Missing columns in data: id.")
})

test_that("convertListCfd2Cfd works", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 0.9), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(1, 5), time = c(0, 0.25, 0.5, 0.9, 1), state = as.factor(c("1_1", "1_2", "2_2", "2_3", "1_3"))
  )
  out <- convertListCfd2Cfd(x)
  expect_equal(out, expectedOut)
})

test_that("convertListCfd2Cfd works with custom sep", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 0.9), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(1, 5), time = c(0, 0.25, 0.5, 0.9, 1), state = as.factor(c("1 & 1", "1 & 2", "2 & 2", "2 & 3", "1 & 3"))
  )
  out <- convertListCfd2Cfd(x, sep = " & ")
  expect_equal(out, expectedOut)
})
