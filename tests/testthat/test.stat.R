# Author: Quentin Grimonprez

context("Statistics on data object")

test_that("compute_time_spent_intern works", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))

  out <- compute_time_spent_intern(dat, 1:3)
  expectedOut <- c(1.5 + 1, 1, 2.5)

  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})


test_that("compute_time_spent_intern works with more labels", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))

  out <- compute_time_spent_intern(dat, 1:4)
  expectedOut <- c(1.5 + 1, 1, 2.5, 0)

  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})

test_that("compute_time_spent works", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))

  out <- compute_time_spent(dat)
  expectedOut <- rbind(
    c(1.5 + 1, 1, 2.5),
    c(3, 3, 0)
  )
  colnames(expectedOut) <- 1:3
  rownames(expectedOut) <- 1:2
  class(expectedOut) <- "timeSpent"

  expect_equal(out, expectedOut)
})

test_that("compute_time_spent keeps unused levels", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  dat$state <- factor(dat$state, levels = 1:4)

  out <- compute_time_spent(dat)
  expectedOut <- rbind(
    c(1.5 + 1, 1, 2.5, 0),
    c(3, 3, 0, 0)
  )
  colnames(expectedOut) <- 1:4
  rownames(expectedOut) <- 1:2
  class(expectedOut) <- "timeSpent"

  expect_equal(out, expectedOut)
})

test_that("boxplot.timeSpent does not produce warnings", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))

  out <- compute_time_spent(dat)

  expect_warning(boxplot(out), regexp = NA)
  expect_warning(boxplot(out, col = c("red", "blue", "green"), outlier.colour = "black"), regexp = NA)
  expect_silent(boxplot(out, col = c("red", "blue", "green"), outlier.colour = "black"))
})


test_that("compute_duration works", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 7, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))

  out <- compute_duration(dat)
  expectedOut <- c("1" = 7, "2" = 6)
  class(expectedOut) <- "duration"

  expect_equivalent(out, expectedOut)
})


test_that("hist.duration does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

  duration <- compute_duration(d_JK)

  expect_warning(hist(duration), regexp = NA)
  expect_warning(hist(duration, color = "red"), regexp = NA)
})


test_that("id_get_state returns the right state", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))

  out <- id_get_state(dat, 3, NAafterTmax = FALSE)
  expect_equal(out, 3)

  out <- id_get_state(dat, 7, NAafterTmax = TRUE)
  expect_equal(out, NA)

  out <- id_get_state(dat, 7, NAafterTmax = FALSE)
  expect_equal(out, 1)

  out <- id_get_state(dat, 6, NAafterTmax = TRUE)
  expect_equal(out, 1)
})


test_that("get_state returns the right state", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))

  expect_error(get_state(dat, c(3, 2)), regexp = "t must be a real.")
  expect_error(get_state(dat, NA), regexp = "t must be a real.")
  expect_error(get_state(dat, NaN), regexp = "t must be a real.")

  out <- get_state(dat, 3)

  expect_equivalent(out, c(3, 2))
})


test_that("estimate_pt works with same t", {
  dat <- data.frame(id = rep(1:2, each = 6), time = rep(0:5, 2), state = c(
    1, 3, 2, 1, 1, 1,
    2, 3, 1, 2, 3, 1
  ))
  out <- estimate_pt(dat)

  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, 0:5)
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1 / 2, 1 / 2, 0, 0, 0, 1, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0, 1 / 2, 0, 1 / 2, 1, 0, 0), nrow = 3, dimnames = list(1:3, 0:5)))
})


test_that("estimate_pt keeps unused levels", {
  dat <- data.frame(id = rep(1:2, each = 6), time = rep(0:5, 2), state = c(
    1, 3, 2, 1, 1, 1,
    2, 3, 1, 2, 3, 1
  ))
  dat$state <- factor(dat$state, levels = 1:4)
  out <- estimate_pt(dat)

  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, 0:5)
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1 / 2, 1 / 2, 0, 0, 0, 0, 1, 0, 1 / 2, 1 / 2, 0, 0, 1 / 2, 1 / 2, 0, 0, 1 / 2, 0, 1 / 2, 0, 1, 0, 0, 0), nrow = 4, dimnames = list(1:4, 0:5)))
})



test_that("estimate_pt works with different t", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 6), state = c(
    1, 3, 2, 1, 1, 1,
    2, 3, 1, 2, 2
  ))
  out <- estimate_pt(dat)

  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, c(0, 1, 1.5, 2, 3, 3.5, 4, 5, 6))
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1 / 2, 1 / 2, 0, 0, 1 / 2, 1 / 2, 0, 0, 1, 1 / 2, 1 / 2, 0, 1, 0, 0, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0), nrow = 3, dimnames = list(1:3, out$t)))


  out <- estimate_pt(dat, NAafterTmax = TRUE)

  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, c(0, 1, 1.5, 2, 3, 3.5, 4, 5, 6))
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1 / 2, 1 / 2, 0, 0, 1 / 2, 1 / 2, 0, 0, 1, 1 / 2, 1 / 2, 0, 1, 0, 0, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2, 0, 0, 1, 0), nrow = 3, dimnames = list(1:3, out$t)))
})


test_that("get_proba returns the right probabilities", {
  pt <- list(t = 1:10, pt = matrix(1:30, nrow = 3, ncol = 10, byrow = TRUE, dimnames = list(letters[1:3], 1:10)))
  class(pt) <- "pt"


  expect_equal(get_proba(pt, 1.5), c("a" = 1, "b" = 11, "c" = 21))
  expect_equal(get_proba(pt, 1), c("a" = 1, "b" = 11, "c" = 21))
  expect_equal(get_proba(pt, 5), c("a" = 5, "b" = 15, "c" = 25))
  expect_equal(get_proba(pt, 5.5), c("a" = 5, "b" = 15, "c" = 25))
  expect_equal(get_proba(pt, 6), c("a" = 6, "b" = 16, "c" = 26))
  expect_equal(get_proba(pt, 10), c("a" = 10, "b" = 20, "c" = 30))
  expect_equal(get_proba(pt, 11), c("a" = 10, "b" = 20, "c" = 30))
  expect_equal(get_proba(pt, 0), c("a" = NA, "b" = NA, "c" = NA))
})


test_that("plot_pt_classic does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

  d_JK2 <- cut_data(d_JK, 10)

  pt <- estimate_pt(d_JK2)

  expect_warning(plot_pt_classic(pt), regexp = NA)
})


test_that("plot_pt_ribbon does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

  d_JK2 <- cut_data(d_JK, 10)

  pt <- estimate_pt(d_JK2)

  expect_warning(plot_pt_ribbon(pt), regexp = NA)
})

test_that("plot_pt does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

  d_JK2 <- cut_data(d_JK, 10)

  pt <- estimate_pt(d_JK2)

  expect_warning(plot(pt, ribbon = FALSE), regexp = NA)
  expect_warning(plot(pt, ribbon = FALSE, col = c("red", "blue", "green", "black")), regexp = NA)
  expect_silent(plot(pt, ribbon = FALSE, col = c("red", "blue", "green", "black")))
  expect_silent(plot(pt, ribbon = TRUE, col = c("red", "blue", "green", "black")))
  expect_warning(plot(pt, ribbon = TRUE, addBorder = TRUE), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE, addBorder = FALSE), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE, addBorder = FALSE, col = c("red", "blue", "green", "black")), regexp = NA)
})


test_that("compute_number_jumps works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- compute_number_jumps(dat)
  expectedOut <- c(5, 4)
  class(expectedOut) <- "njump"

  expect_equivalent(out, expectedOut)
})

test_that("compute_number_jumpsIntern works with countDuplicated = TRUE", {
  dat <- data.frame(id = 1:20, time = 1:20, state = rep(1:10, each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = TRUE)
  expectedOut <- 19

  expect_equivalent(out, expectedOut)


  dat <- data.frame(id = 1:20, time = 1:20, state = rep(letters[1:10], each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = TRUE)
  expectedOut <- 19

  expect_equivalent(out, expectedOut)
})

test_that("compute_number_jumpsIntern works with countDuplicated = FALSE", {
  # without duplicate state
  dat <- data.frame(id = 1:20, time = 1:20, state = 1:20)
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE)
  expectedOut <- 19

  expect_equivalent(out, expectedOut)

  # with ordered time
  dat <- data.frame(id = 1:20, time = 1:20, state = rep(1:10, each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE)
  expectedOut <- 9

  expect_equivalent(out, expectedOut)

  # with unordered time
  dat <- data.frame(id = 1:20, time = c(11:20, 1:10), state = as.factor(rep(letters[1:5], each = 4)))
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE)
  expectedOut <- 5

  expect_equivalent(out, expectedOut)
})



test_that("hist.njump does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

  njump <- compute_number_jumps(d_JK)

  expect_warning(hist(njump), regexp = NA)
  expect_warning(hist(njump, color = "red"), regexp = NA)
})


test_that("statetable works", {
  dat <- data.frame(id = rep(1:2, c(6, 6)), time = c(0:5, 0, 1.5, 2, 3.5, 5, 6), state = c(1:6, 1:5, 5))
  out <- statetable(dat, removeDiagonal = FALSE)

  expectedOut <- matrix(0, nrow = 6, ncol = 6)
  expectedOut[1, 2] <- expectedOut[2, 3] <- expectedOut[3, 4] <- expectedOut[4, 5] <- 2
  expectedOut[5, 6] <- expectedOut[5, 5] <- 1

  expect_equivalent(out, expectedOut)


  out <- statetable(dat, removeDiagonal = TRUE)
  expectedOut[5, 5] <- 0

  expect_equivalent(out, expectedOut)
})


test_that("rep_large_ind works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- rep_large_ind(dat)
  expectedOut <- data.frame(
    id = rep(1:2, c(5, 4)),
    t_start = c(0:4, 0, 1.5, 2, 3.5),
    t_end = c(1:5, 1.5, 2, 3.5, 5),
    state = c(1:5, 1:4)
  )

  expect_equivalent(out, expectedOut)



  out <- rep_large_ind(dat[dat$id == 1, ])

  expect_equivalent(out, expectedOut[expectedOut$id == 1, ])
})


test_that("orderFirstState works", {
  dat <- data.frame(id = rep(1:5, c(2, 1, 2, 2, 1)), time = c(0:1, 0, 0, 2, 0, 3, 0), state = c(1:2, 1, 1:2, 2:1, 2))
  expectedOut <- data.frame(id = c(1, 3, 2, 4, 5), time = c(1, 2, Inf, 3, Inf), state = c(1, 1, 1, 2, 2))

  out <- orderFirstState(dat)

  expect_equivalent(out, expectedOut)
})


test_that("computePosition works", {
  dat <- data.frame(id = rep(1:5, each = 2), time = c(0:1, c(0, 2), c(0, 1), c(0, 3), c(0, 2)), state = c(1:2, 2:1, 2:1, 2:1, 1:2))
  d <- rep_large_ind(dat)

  out <- computePosition(dat, d$id, sort = FALSE)

  expect_equivalent(out, d$id)


  out <- computePosition(dat, d$id, sort = TRUE)

  expect_equivalent(out, order(c(1, 5, 3, 2, 4)))
})


test_that("computePositionPerGroup works", {
  dat <- data.frame(
    id = rep(1:5, each = 2), time = c(0:1, c(0, 2), c(0, 1), c(0, 3), c(0, 2)),
    state = c(1:2, 2:1, 2:1, 2:1, 1:2), group = rep(1:2, c(6, 4))
  )
  d <- rep_large_ind(dat)

  out <- computePositionPerGroup(dat, d$id, d$group, sort = FALSE)

  expect_equivalent(out, d$id)

  # non consecutive group number
  d$group[d$group == 2] <- 3
  dat$group[dat$group == 2] <- 3
  out <- computePositionPerGroup(dat, d$id, d$group, sort = TRUE)

  expect_equivalent(out, c(1, 3, 2, 5, 4))
})

test_that("computePositionPerGroup works when there is only 1 group", {
  dat <- data.frame(
    id = rep(1:5, each = 2), time = c(0:1, c(0, 2), c(0, 1), c(0, 3), c(0, 2)),
    state = c(1:2, 2:1, 2:1, 2:1, 1:2), group = rep(1, 10)
  )
  d <- rep_large_ind(dat)

  out <- computePositionPerGroup(dat, d$id, d$group, sort = FALSE)

  expect_equivalent(out, d$id)
})


test_that("createLabeller works", {
  group <- rep(1:3, 6:4)
  f <- createLabeller(group)

  expect_is(f, "function")
  expect_equal(f(value = "1"), list("1" = "1: n=6"))
  expect_equal(f(value = "2"), list("2" = "2: n=5"))
  expect_equal(f(value = "3"), list("3" = "3: n=4"))
})


test_that("plotData does not produce warnings", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
  d_JKT <- cut_data(d_JK, Tmax = 10)
  group <- rep(1:2, c(3, 7))

  expect_warning(plotData(d_JK, addId = TRUE, addBorder = TRUE, sort = FALSE), regexp = NA)
  expect_warning(plotData(d_JK, addId = FALSE, addBorder = FALSE, col = c("red", "blue", "green", "yellow")), regexp = NA)
  expect_silent(plotData(d_JK, addId = FALSE, addBorder = FALSE, col = c("red", "blue", "green", "yellow")))
  expect_warning(plotData(d_JK, addId = FALSE, addBorder = FALSE, sort = TRUE), regexp = NA)
  expect_warning(plotData(d_JK, group = group, addId = FALSE, addBorder = FALSE, sort = FALSE), regexp = NA)
  expect_warning(plotData(d_JK, group = group, addId = FALSE, addBorder = FALSE, sort = TRUE), regexp = NA)
  expect_warning(plotData(d_JK, group = group, addId = FALSE, addBorder = FALSE, sort = TRUE, nCol = 2), regexp = NA)
  expect_warning(plotData(d_JK, group = as.factor(group), addId = FALSE, addBorder = FALSE, sort = TRUE, nCol = 2), regexp = NA)
})

test_that("plotData produces an error when group is bad", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
  d_JKT <- cut_data(d_JK, Tmax = 10)


  expect_error(plotData(d_JK, group = 2, addId = FALSE, addBorder = FALSE, sort = FALSE), regexp = "group must be a vector with the same length than the number of ids of data.")
  expect_error(plotData(d_JK, group = 2:nrow(d_JK), addId = FALSE, addBorder = FALSE, sort = FALSE), regexp = "group must be a vector with the same length than the number of ids of data.")
})

test_that("plotData produces an error when nCol is bad", {
  K <- 4
  PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
  lambda_PJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
  d_JKT <- cut_data(d_JK, Tmax = 10)


  expect_error(plotData(d_JK, group = rep(1:2, each = 5), addId = FALSE, addBorder = FALSE, sort = FALSE, nCol = -1), regexp = "nCol must be an integer > 0.")
  expect_error(plotData(d_JK, group = rep(1:2, each = 5), addId = FALSE, addBorder = FALSE, sort = FALSE, nCol = "aaa"), regexp = "nCol must be an integer > 0.")
  expect_error(plotData(d_JK, group = rep(1:2, each = 5), addId = FALSE, addBorder = FALSE, sort = FALSE, nCol = 1:3), regexp = "nCol must be an integer > 0.")
})

test_that("summary_cfd words", {
  dat <- data.frame(id = rep(1:5, c(2, 1, 2, 2, 1)), time = c(0:1, 0, 0, 2, 0, 3, 0), state = c(1:2, 1, 1:2, 2:1, 2))

  expect_output(out <- summary_cfd(dat))

  expectedOut <- list(
    nRow = 8, nInd = 5, timeRange = c(0, 3), uniqueStart = TRUE, uniqueEnd = FALSE,
    states = c("1", "2"), visit = array(c(4, 4), dimnames = list(c("1", "2")))
  )

  expect_equal(out, expectedOut)
})
