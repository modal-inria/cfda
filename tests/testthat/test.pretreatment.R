# Author: Quentin Grimonprez

context("Data pretreatment")

test_that("cut_cfd with equal Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 4, prolongLastState = "all")
  expect_equal(out, dat)

  out <- cut_cfd(dat, Tmax = 4, prolongLastState = c())
  expect_equal(out, dat)
})

test_that("cut_cfd with lower Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 3, prolongLastState = "all")
  expectedOut <- dat
  expectedOut[3, 2:3] <- c(3, 3)

  expect_equal(out, expectedOut)

  out <- cut_cfd(dat, Tmax = 3, prolongLastState = c())

  expect_equal(out, expectedOut)
})

test_that("cut_cfd with lower Tmax and a time value equal to the desired Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 2, prolongLastState = "all")
  expectedOut <- dat[1:2, ]

  expect_equal(out, expectedOut)

  out <- cut_cfd(dat, Tmax = 2, prolongLastState = c())
  expect_equal(out, expectedOut)
})

test_that("cut_cfd with greater Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 5, prolongLastState = "all")
  expectedOut <- dat
  expectedOut[4, 1:3] <- c(1, 5, 2)

  expect_equal(out, expectedOut)

  expect_warning(out <- cut_cfd(dat, Tmax = 5, prolongLastState = c()),
                 regexp = paste0("id ", 1, " does not end with an absorbing state. Cannot impute the state until time ",
                                 5, ". Please, add more records or change the Tmax value."))
  expectedOut$state[3:4] <- NA
  expect_equal(out, expectedOut)
})


test_that("cut_data generates error with bad parameters", {
  dat <- data.frame(id = rep(1:3, each = 3), time = c(0, 2, 4, 0, 1.5, 5, 0, 2.5, 3), state = c(1, 3, 2, 1, 2, 3, 1, 3, 1))

  expect_error(cut_data(dat, Tmax = c(4, 5), prolongLastState = "all"), regexp = "Tmax must be a real.")
  expect_error(cut_data(dat, Tmax = NA, prolongLastState = "all"), regexp = "Tmax must be a real.")
  expect_error(cut_data(dat, Tmax = NaN, prolongLastState = "all"), regexp = "Tmax must be a real.")
  expect_error(cut_data(dat, Tmax = 2, prolongLastState = "all", warning = "fff"),
               regexp = "warning must be either TRUE or FALSE.")
  expect_error(cut_data(dat, Tmax = 2, prolongLastState = "all", warning = TRUE, NAstate = c(1, 2)),
               regexp = "NAstate must have a length of 1")
})

test_that("cut_data works", {
    dat <- data.frame(id = rep(1:3, each = 3), time = c(0, 2, 4, 0, 1.5, 5, 0, 2.5, 3), state = c(1, 3, 2, 1, 2, 3, 1, 3, 1))

  out <- cut_data(dat, Tmax = 4, prolongLastState = "all")
  expectedOut <- dat
  expectedOut[6, 1:3] <- c(2, 4, 2)
  expectedOut[10, 1:3] <- c(3, 4, 1)

  expect_equivalent(out, expectedOut)

  expect_warning(out <- cut_data(dat, Tmax = 4, prolongLastState = c(), warning = TRUE),
                 regexp = paste0("id ", 3, " does not end with an absorbing state. Cannot impute the state until time ",
                                 4, ". Please, add more records or change the Tmax value."))
  expectedOut$state[9:10] <- "Not observable"
  expect_equal(out, expectedOut)

  expect_silent(out <- cut_data(dat, Tmax = 4, prolongLastState = c(), NAstate = "blabla"))
  expectedOut$state[9:10] <- "blabla"
  expect_equal(out, expectedOut)
})


test_that("refactorCategorical works when oldCateg and newCateg do not have common elements", {
  x <- letters[c(26:9, 1:8, 20:25)]
  oldCateg <- letters
  newCateg <- seq_along(oldCateg)

  expectedOut <- c(26:9, 1:8, 20:25)

  out <- refactorCategorical(x, oldCateg, newCateg)
  expect_equal(as.character(out), as.character(expectedOut))
})

test_that("refactorCategorical works when oldCateg and newCateg have common elements", {
  x <- as.character(c(7:10, 0:6, 7:10))
  oldCateg <- as.character(0:10)
  newCateg <- 1:11

  expectedOut <- c(7:10, 0:6, 7:10) + 1

  out <- refactorCategorical(x, oldCateg, newCateg)
  expect_equal(as.character(out), as.character(expectedOut))
})

test_that("refactorCategorical works when some categories are merged", {
  x <- letters[c(1, 4:6, 2:3)]
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))
  expectedOut <- c("voyelle", "consonne", "voyelle", "consonne", "consonne", "consonne")

  out <- refactorCategorical(x, oldCateg, newCateg)

  expect_equal(as.character(out), as.character(expectedOut))
})

test_that("refactorCategorical works when there are categories not included in the data", {
  x <- letters[c(1, 4:6, 2:3)]
  oldCateg <- letters[1:7]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 2))

  expectedOut <- c("voyelle", "consonne", "voyelle", "consonne", "consonne", "consonne")

  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg), regexp = NA)
  expect_equal(as.character(out), as.character(expectedOut))


  x <- letters[c(1, 4:7, 2:3)]
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))

  expectedOut <- c("voyelle", "consonne", "voyelle", "consonne", NA, "consonne", "consonne")

  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg))
  expect_equal(as.character(out), as.character(expectedOut))
})


test_that("refactorCategorical kept NA values in data", {
  x <- c(letters[c(1, 4:6, 2:3)], NA)
  oldCateg <- letters[1:6]
  newCateg <- rep(c("voyelle", "consonne", "voyelle", "consonne"), c(1, 3, 1, 1))

  expectedOut <- c("voyelle", "consonne", "voyelle", "consonne", "consonne", "consonne", NA)

  expect_warning(out <- refactorCategorical(x, oldCateg, newCateg), regexp = NA)
  expect_equal(as.character(out), as.character(expectedOut))
})


test_that("stateToInteger works", {
  x <- letters[c(1, 4:6, 2:3)]

  out <- stateToInteger(x)
  expectedOut <- list(state = c(1, 4:6, 2:3), label = data.frame(label = sort(x), code = 1:6))

  expect_equal(out, expectedOut)


  x <- letters[c(6:1, 6, 2)]

  out <- stateToInteger(x)
  expectedOut <- list(state = c(6:1, 6, 2), label = data.frame(label = letters[1:6], code = 1:6))

  expect_equal(out, expectedOut)
})


test_that("remove_duplicated_states.intern works with keep.last = FALSE", {
  data <- data.frame(id = rep(1, 10), time = 1:10, state = rep(1:5, each = 2))

  out <- remove_duplicated_states.intern(data, keep.last = FALSE)
  expectedOut <- data.frame(id = rep(1, 5), time = 1:5 * 2 - 1, state = 1:5)

  expect_equivalent(out, expectedOut)


  data$state <- as.factor(data$state)

  out <- remove_duplicated_states.intern(data, keep.last = FALSE)
  expectedOut <- data.frame(id = rep(1, 5), time = 1:5 * 2 - 1, state = as.factor(1:5))

  expect_equivalent(out, expectedOut)
})


test_that("remove_duplicated_states.intern works with keep.last = TRUE", {
  data <- data.frame(id = rep(1, 10), time = 1:10, state = rep(1:5, each = 2))

  out <- remove_duplicated_states.intern(data, keep.last = TRUE)
  expectedOut <- data.frame(id = rep(1, 6), time = c(1:5 * 2 - 1, 10), state = c(1:5, 5))

  expect_equivalent(out, expectedOut)
})

test_that("remove_duplicated_states works", {
  data <- data.frame(id = rep(1:3, c(10, 3, 8)), time = c(1:10, 1:3, 1:8), state = c(rep(1:5, each = 2), 1:3, rep(1:3, c(1, 6, 1))))

  out <- remove_duplicated_states(data, keep.last = FALSE)
  expectedOut <- data.frame(id = rep(1:3, c(5, 3, 3)), time = c(1:5 * 2 - 1, 1:3, 1, 2, 8), state = c(1:5, 1:3, 1:3))

  expect_equivalent(out, expectedOut)
})

test_that("matrixToCfd works", {
  x <- matrix(c("a", "b", "c", "c",
                "c", "a", "a", "a",
                "b", "c", "a", "b"), ncol = 4, byrow = TRUE)


  out <- matrixToCfd(x, byrow = FALSE)

  expectedOut <- data.frame(id = rep(1:4, each = 3),
                            time = rep(1:3, 4),
                            state = c("a", "c", "b", "b", "a", "c", "c", "a", "a", "c", "a", "b"))

  expect_equivalent(out, expectedOut)


  out <- matrixToCfd(x, byrow = TRUE)

  expectedOut <- data.frame(id = rep(1:3, c(4, 3, 4)),
                            time = c(1:4, c(1, 2, 4), 1:4),
                            state = c("a", "b", "c", "c",
                                      "c", "a", "a",
                                      "b", "c", "a", "b"))

  expect_equivalent(out, expectedOut)

  out <- matrixToCfd(x, times = c(1.5, 2.5, 3.5), byrow = FALSE)

  expectedOut <- data.frame(id = rep(1:4, each = 3),
                            time = rep(1:3, 4) + 0.5,
                            state = c("a", "c", "b", "b", "a", "c", "c", "a", "a", "c", "a", "b"))

  expect_equivalent(out, expectedOut)


  out <- matrixToCfd(x, times = c(1.5, 2.5, 3.5, 4.5), byrow = TRUE)

  expectedOut <- data.frame(id = rep(1:3, c(4, 3, 4)),
                            time = c(1:4, c(1, 2, 4), 1:4) + 0.5,
                            state = c("a", "b", "c", "c",
                                      "c", "a", "a",
                                      "b", "c", "a", "b"))

  expect_equivalent(out, expectedOut)

  times <- matrix(c(1:3, 0:2 + 0.5, 1:3, 1:3 + 0.5), nrow = 3)
  out <- matrixToCfd(x, times = times, byrow = FALSE)
  expectedOut <- data.frame(id = rep(1:4, each = 3),
                            time = as.vector(times),
                            state = c("a", "c", "b", "b", "a", "c", "c", "a", "a", "c", "a", "b"))
  expect_equivalent(out, expectedOut)

})

test_that("matrixToCfd keeps manages labels", {
  x <- matrix(c("a", "b", "c", "c",
                "c", "a", "a", "a",
                "b", "c", "a", "b"), ncol = 4, byrow = TRUE,
              dimnames = list(paste0("time", 1:3), paste0("ind", 1:4)))


  out <- matrixToCfd(x, byrow = FALSE)

  expectedOut <- data.frame(id = rep(paste0("ind", 1:4), each = 3),
                            time = rep(1:3, 4),
                            state = c("a", "c", "b", "b", "a", "c", "c", "a", "a", "c", "a", "b"))

  expect_equal(out, expectedOut)

  x <- matrix(c("a", "b", "c", "c",
                "c", "a", "a", "a",
                "b", "c", "a", "b"), ncol = 4, byrow = TRUE,
              dimnames = list(paste0("time", 1:3), paste0("ind", 1:4)))


  out <- matrixToCfd(x, byrow = FALSE, labels = c("a", "b", "c", "d"))

  expectedOut <- data.frame(id = rep(c("a", "b", "c", "d"), each = 3),
                            time = rep(1:3, 4),
                            state = c("a", "c", "b", "b", "a", "c", "c", "a", "a", "c", "a", "b"))

  expect_equal(out, expectedOut)
})


test_that("matrixToCfd errors", {
  x <- matrix(c("a", "b", "c", "c",
                "c", "a", "a", "a",
                "b", "c", "a", "b"), ncol = 4, byrow = TRUE)

  expect_error(matrixToCfd(x, times = NULL, byrow = 3), "byrow must be either TRUE or FALSE.")
  expect_error(matrixToCfd(x, times = 3, byrow = TRUE), "times must be a numeric vector of length 4")
  expect_error(matrixToCfd(x, times = c(3, 2), byrow = FALSE), "times must be a numeric vector of length 3")
  expect_error(matrixToCfd(x, times = c("a", "b", "c"), byrow = FALSE), "times must be a numeric vector of length 3")
  expect_error(matrixToCfd(c(1, 3), times = NULL, byrow = TRUE), "X must be a matrix or a data.frame")
  expect_error(matrixToCfd(x, times = NULL, labels = c("a"), byrow = TRUE), "labels must be a vector of length 3")
})


data("CanadianWeather")
temp <- CanadianWeather$dailyAv[,, "Temperature.C"]
basis <- create.bspline.basis(c(1, 365), nbasis = 8, norder = 4)
fd <- smooth.basis(1:365, temp, basis)$fd


test_that("quantiMatrixToCfd works", {
  out <- quantiMatrixToCfd(temp, breaks = c(-50, -10, 0, 10, 20, 50), right = FALSE,
                          labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), idLabels = NULL, times = 0:364)

  expect_true(is.data.frame(out))
  expect_equal(colnames(out), c("id", "time", "state"))
  expect_equal(range(out$time), c(0, 364))
  expect_equal(sort(unique(out$state)), sort(c("Very Cold", "Cold", "Fresh", "OK", "Hot")))
  expect_equal(unique(out$id), colnames(temp))

  expectedOut <- data.frame(id = rep("St. Johns", 10),
                            time = c(0, 94, 164, 271, 275, 276, 335, 340, 341, 364),
                            state = c("Cold", "Fresh", "OK", "Fresh", "OK", "Fresh", "Cold", "Fresh", "Cold", "Cold"))
  expect_equal(out[out$id==out$id[1], ], expectedOut)
})

test_that("quantiMatrixToCfd errors", {
  expect_error(quantiMatrixToCfd(temp, breaks = c(-50, -10), right = TRUE,
                                labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365))
  expect_error(quantiMatrixToCfd(temp, breaks = c(-50, -10), right = "3",
                                labels = c("Very Cold"), times = 1:365))
  expect_error(quantiMatrixToCfd(fd, breaks = c(-50, -10), right = TRUE,
                                labels = c("Very Cold"), times = 1:365))
})

test_that("fdToCfd works", {
  out <- fdToCfd(fd, breaks = c(-50, -10, 0, 10, 20, 50), right = FALSE,
                 labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365)

  expect_true(is.data.frame(out))
  expect_equal(colnames(out), c("id", "time", "state"))
  expect_equal(range(out$time), c(1, 365))
  expect_equal(sort(unique(out$state)), sort(c("Very Cold", "Cold", "Fresh", "OK", "Hot")))
  expect_equal(unique(out$id), colnames(temp))

  expectedOut <- data.frame(id = rep("St. Johns", 6),
                            time = c(1, 97, 161, 271, 340, 365),
                            state = c("Cold", "Fresh", "OK", "Fresh", "Cold", "Cold"))
  expect_equal(out[out$id==out$id[1], ], expectedOut)
})

test_that("fdToCfd errors", {
  expect_error(fdToCfd(fd, breaks = c(-50, -10), right = TRUE,
                       labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365))
  expect_error(fdToCfd(fd, breaks = c(-50, -10), right = "3",
                       labels = c("Very Cold"), times = 1:365))
  expect_error(fdToCfd(5, breaks = c(-50, -10), right = TRUE,
                       labels = c("Very Cold"), times = 1:365))
})


test_that("convertToCfd works with matrix", {
  out <- convertToCfd(temp, breaks = c(-50, -10, 0, 10, 20, 50), right = FALSE,
                      labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), idLabels = NULL, times = 0:364)

  expect_true(is.data.frame(out))
  expect_equal(colnames(out), c("id", "time", "state"))
  expect_equal(range(out$time), c(0, 364))
  expect_equal(sort(unique(out$state)), sort(c("Very Cold", "Cold", "Fresh", "OK", "Hot")))
  expect_equal(unique(out$id), colnames(temp))

  expectedOut <- data.frame(id = rep("St. Johns", 10),
                            time = c(0, 94, 164, 271, 275, 276, 335, 340, 341, 364),
                            state = c("Cold", "Fresh", "OK", "Fresh", "OK", "Fresh", "Cold", "Fresh", "Cold", "Cold"))
  expect_equal(out[out$id==out$id[1], ], expectedOut)
})

test_that("convertToCfd works with fd", {
  out <- convertToCfd(fd, breaks = c(-50, -10, 0, 10, 20, 50), right = FALSE,
                      labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365)

  expect_true(is.data.frame(out))
  expect_equal(colnames(out), c("id", "time", "state"))
  expect_equal(range(out$time), c(1, 365))
  expect_equal(sort(unique(out$state)), sort(c("Very Cold", "Cold", "Fresh", "OK", "Hot")))
  expect_equal(unique(out$id), colnames(temp))

  expectedOut <- data.frame(id = rep("St. Johns", 6),
                            time = c(1, 97, 161, 271, 340, 365),
                            state = c("Cold", "Fresh", "OK", "Fresh", "Cold", "Cold"))
  expect_equal(out[out$id==out$id[1], ], expectedOut)
})

test_that("convertToCfd errors", {
  expect_error(convertToCfd(fd, breaks = c(-50, -10), right = TRUE,
                            labels = c("Very Cold"), times = 1:365),
               "The conversion has generated NA. Please, correct your breaks.")

  expect_error(convertToCfd(temp, breaks = c(-50, -10), right = TRUE,
                            labels = c("Very Cold"), times = 1:365),
               "The conversion has generated NA. Please, correct your breaks.")
})

