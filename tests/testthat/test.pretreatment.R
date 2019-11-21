# Author: Quentin Grimonprez

context("Data pretreatment")

test_that("cut_cfd with equal Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))
  
  out <- cut_cfd(dat, Tmax = 4)
  
  expect_equal(out, dat)
})

test_that("cut_cfd with lower Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))
  
  out <- cut_cfd(dat, Tmax = 3)
  expectedOut <- dat
  expectedOut[3, 2:3] = c(3, 3)
  
  expect_equal(out, expectedOut)
})

test_that("cut_cfd with greater Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))
  
  out <- cut_cfd(dat, Tmax = 5)
  expectedOut <- dat
  expectedOut[4, 1:3] = c(1, 5, 2)
  
  expect_equal(out, expectedOut)
})


test_that("msm2msmTmax works", {
  dat <- data.frame(id = rep(1:3, each = 3), time = c(0, 2, 4, 0, 1.5, 5, 0, 2.5, 3), state = c(1, 3, 2, 1, 2, 3, 1, 3, 1))
  
  out <- msm2msmTmax(dat, Tmax = 4)
  expectedOut <- dat
  expectedOut[6, 1:3] = c(2, 4, 2)
  expectedOut[10, 1:3] = c(3, 4, 1)
  
  expect_equivalent(out, expectedOut)
})


test_that("refactorCategorical works when oldCateg and newCateg do not have common elements", {
  x <- letters
  oldCateg <- letters
  newCateg <- seq_along(oldCateg)
  
  expectedOut <- newCateg
  
  out <- refactorCategorical(x,oldCateg,newCateg)
  expect_equal(out, expectedOut)
})

test_that("refactorCategorical works when oldCateg and newCateg have common elements", {
  x <- as.character(0:10)
  oldCateg <- as.character(0:10)
  newCateg <- 1:11
  
  expectedOut <- newCateg
  
  out <- refactorCategorical(x, oldCateg, newCateg)
  expect_equal(out, expectedOut)
})

test_that("refactorCategorical works when some categories are merged", {
  
  x <- letters[1:6]
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))
  expectedOut <- newCateg
  
  out <- refactorCategorical(x, oldCateg, newCateg)
  
  expect_equal(out, expectedOut)
})

test_that("refactorCategorical works when there are categories not included in the data", {
  
  x <- letters[1:6]
  oldCateg <- letters[1:7]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 2))
  
  expectedOut <- newCateg[1:6]
  
  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg), regexp = NA)
  expect_equal(out, expectedOut)
  
  
  x <- letters[1:7]
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))
  
  expectedOut <- c(newCateg[1:6], NA)
  
  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg))
  expect_equal(out, expectedOut)
})


test_that("refactorCategorical kept NA values in data", {
  
  x <- c(letters[1:6], NA)
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))
  expectedOut <- c(newCateg[1:6], NA)
  
  
  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg), regexp = NA)
  expect_equal(out, expectedOut)
})


test_that("stateToInteger works", {
  x <- letters[1:6]

  out <- stateToInteger(x)
  expectedOut <- list(state = 1:6, label = data.frame(label = x, code = 1:6))
  
  expect_equal(out, expectedOut)
  
  
  x <- letters[6:1]
  
  out <- stateToInteger(x)
  expectedOut <- list(state = 6:1, label = data.frame(label = letters[1:6], code = 1:6))
  
  expect_equal(out, expectedOut)
})