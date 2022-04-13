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
  expectedOut[3, 2:3] <- c(3, 3)

  expect_equal(out, expectedOut)
})

test_that("cut_cfd with lower Tmax and a time value equal to the desired Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 2)
  expectedOut <- dat[1:2, ]

  expect_equal(out, expectedOut)
})

test_that("cut_cfd with greater Tmax", {
  dat <- data.frame(id = rep(1, 3), time = c(0, 2, 4), state = c(1, 3, 2))

  out <- cut_cfd(dat, Tmax = 5)
  expectedOut <- dat
  expectedOut[4, 1:3] <- c(1, 5, 2)

  expect_equal(out, expectedOut)
})


test_that("cut_data works", {
  dat <- data.frame(id = rep(1:3, each = 3), time = c(0, 2, 4, 0, 1.5, 5, 0, 2.5, 3), state = c(1, 3, 2, 1, 2, 3, 1, 3, 1))

  expect_error(cut_data(dat, Tmax = c(4, 5)), regexp = "Tmax must be a real.")
  expect_error(cut_data(dat, Tmax = NA), regexp = "Tmax must be a real.")
  expect_error(cut_data(dat, Tmax = NaN), regexp = "Tmax must be a real.")

  out <- cut_data(dat, Tmax = 4)
  expectedOut <- dat
  expectedOut[6, 1:3] <- c(2, 4, 2)
  expectedOut[10, 1:3] <- c(3, 4, 1)

  expect_equivalent(out, expectedOut)
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


test_that("quanti2quali works", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE)

  out <- quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = TRUE, labels = NULL)
  expectedOut <- matrix(c(2, 2, 3, 4,
                          2, 3, 3, 3,
                          2, 3, 4, 3), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)

  out <- quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = FALSE, labels = NULL)
  expectedOut <- matrix(c(2, 2, 3, 3,
                          1, 2, 3, 3,
                          2, 3, 4, 3), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)
})

test_that("quanti2quali manages the first/last interval as closed", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE)

  out <- quanti2quali(x, c(0, 1, 2.2), leftClosed = TRUE, labels = NULL)
  expectedOut <- matrix(c(1, 1, 2, 2,
                          1, 2, 2, 2,
                          1, 2, 2, 2), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)

  out <- quanti2quali(x, c(0, 1, 2.2), leftClosed = FALSE, labels = NULL)
  expectedOut <- matrix(c(1, 1, 2, 2,
                          1, 1, 2, 2,
                          1, 2, 2, 2), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)
})

test_that("quanti2quali state labels", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE)

  out <- quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = TRUE, labels = c("a", "b", "c", "d"))
  expectedOut <- matrix(c("b", "b", "c", "d",
                          "b", "c", "c", "c",
                          "b", "c", "d", "c"), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)

  out <- quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = FALSE, labels = c(4, 3, 2, 1))
  expectedOut <- matrix(c(3, 3, 2, 2,
                          4, 3, 2, 2,
                          3, 2, 1, 2), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)
})

test_that("quanti2quali keeps column and row names", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE,
              dimnames = list(paste0("time", 1:3), paste0("ind", 1:4)))

  out <- quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = TRUE, labels = c("a", "b", "c", "d"))
  expectedOut <- matrix(c("b", "b", "c", "d",
                          "b", "c", "c", "c",
                          "b", "c", "d", "c"), ncol = 4, byrow = TRUE,
                        dimnames = list(paste0("time", 1:3), paste0("ind", 1:4)))
  expect_equal(out, expectedOut)
})

test_that("quanti2quali warnings NA elements", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE)

  expect_warning(out <- quanti2quali(x, c(0.2, 1, 2, Inf), leftClosed = FALSE, labels = NULL),
                 "The conversion has generated NA elements")
  expectedOut <- matrix(c(1, 1, 2, 2,
                          NA, 1, 2, 2,
                          1, 2, 3, 2), ncol = 4, byrow = TRUE)
  expect_equivalent(out, expectedOut)
})


test_that("quanti2quali errors", {
  x <- matrix(c(0.5, 0.7, 1.5, 2,
                0, 1, 1.2, 1.5,
                0.8, 1.6, 2.2, 1.7), ncol = 4, byrow = TRUE)

  expect_error(quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = 3, labels = NULL),
               "leftClosed must be either TRUE or FALSE.")
  expect_error(quanti2quali(4, c(-Inf, 0, 1, 2, Inf), leftClosed = TRUE, labels = NULL), "X must be a matrix")
  expect_error(quanti2quali(x, c(-Inf, 0, 1, 2, Inf), leftClosed = TRUE, labels = c(1)), "labels must be a vector of length 4")
})


data("CanadianWeather")
temp <- CanadianWeather$dailyAv[,, "Temperature.C"]
basis <- create.bspline.basis(c(1, 365), nbasis = 8, norder = 4)
fd <- smooth.basis(1:365, temp, basis)$fd


test_that("qualiMatrixToCfd works", {
  out <- qualiMatrixToCfd(temp, thr = c(-50, -10, 0, 10, 20, 50), leftClosed = TRUE,
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

test_that("qualiMatrixToCfd errors", {
  expect_error(qualiMatrixToCfd(temp, thr = c(-50, -10), leftClosed = TRUE,
                                labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365),
               "labels must be a vector of length 1")
  expect_error(qualiMatrixToCfd(temp, thr = c(-50, -10), leftClosed = 3,
                                labels = c("Very Cold"), times = 1:365),
               "leftClosed must be either TRUE or FALSE.")
  expect_error(qualiMatrixToCfd(fd, thr = c(-50, -10), leftClosed = TRUE,
                                labels = c("Very Cold"), times = 1:365),
               "X must be a matrix")
})

test_that("fdToCfd works", {
  out <- fdToCfd(fd, thr = c(-50, -10, 0, 10, 20, 50), leftClosed = TRUE,
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
  expect_error(fdToCfd(fd, thr = c(-50, -10), leftClosed = TRUE,
                       labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"), times = 1:365),
               "labels must be a vector of length 1")
  expect_error(fdToCfd(fd, thr = c(-50, -10), leftClosed = 3,
                       labels = c("Very Cold"), times = 1:365),
               "leftClosed must be either TRUE or FALSE.")
  expect_error(fdToCfd(5, thr = c(-50, -10), leftClosed = TRUE,
                       labels = c("Very Cold"), times = 1:365),
               "fd is not a fd object")
})


test_that("convertToCfd works with matrix", {
  out <- convertToCfd(temp, thr = c(-50, -10, 0, 10, 20, 50), leftClosed = TRUE,
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
  out <- convertToCfd(fd, thr = c(-50, -10, 0, 10, 20, 50), leftClosed = TRUE,
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
