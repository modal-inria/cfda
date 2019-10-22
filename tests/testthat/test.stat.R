context("Statistics on data_msm object")

test_that("compute_Time_Spent_intern works", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_Time_Spent_intern(dat, K = 3)
  expectedOut <- c(1.5 + 1, 1, 2.5)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})


test_that("compute_Time_Spent_intern works with a higher K", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_Time_Spent_intern(dat, K = 4)
  expectedOut <- c(1.5 + 1, 1, 2.5, 0)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})
