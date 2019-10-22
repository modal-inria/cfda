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
