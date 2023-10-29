# Author: Quentin Grimonprez

context("Multivariate Encoding")

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
