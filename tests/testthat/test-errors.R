#
# Test the checks on the parameters when calling caRamel
#

context("errors")

fonseca <- function(i) {
  s2 <- 1 / sqrt(2)
  val1 <- 1 - exp(-(x[i,1] - s2) * (x[i,1] - s2) - (x[i,2] - s2) * (x[i,2] - s2))
  s2 <- 1 / sqrt(2)
  val2 <- 1 - exp(-(x[i,1] + s2) * (x[i,1] + s2) - (x[i,2] + s2) * (x[i,2] + s2))
  return(c(val1, val2))
}

nvar <- 2 # number of variables
bounds <- matrix(data = 1, nrow = nvar, ncol = 2) # upper and lower bounds
bounds[, 1] <- -4 * bounds[, 1]
bounds[, 2] <- 4 * bounds[, 2]
nobj <- 2 # number of objectives
minmax <- c(FALSE, FALSE) # min and min
popsize <- 100 # size of the genetic population
archsize <- 100 # size of the archive for the Pareto front
maxrun <- 1000 # maximum number of calls
prec <- matrix(1.e-3, nrow = 1, ncol = nobj) # accuracy for the convergence phase

test_that("Right number of objectives", {
    # must be greater than one
    nobj <- 1
    results <-
      caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
    nobj <- 2
})

test_that("Right number of variables", {
    # must be strictly positive"
    nvar <- 0
    results <-
      caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
    nvar <- 2
})

test_that("Test the nmuber of goals", {
    # must be equal to nobj
    minmax <- c(TRUE)
    results <- caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
    minmax <- c(FALSE, FALSE)
})

test_that("Test bound arrays", {
    # must be equal to nvar
    fbounds <- matrix(data = 1, nrow = nvar+1, ncol = 2)
    results <- caRamel(nobj,
              nvar,
              minmax,
              fbounds,
              fonseca,
              popsize,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
})

test_that("Test objective function", {
    # must be R function
    results <- caRamel(nobj,
              nvar,
              minmax,
              bounds,
              "kursawe",
              popsize,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
})

test_that("Test size of the genetic population", {
    # must be positive
    results <- caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              0,
              archsize,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
})

test_that("Test size of the archive", {
    # must be positive
    results <- caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              0,
              maxrun,
              prec)
    expect_true(results$success==FALSE)
})

test_that("Test size of the maximum number of runs", {
    # must be strictly positive
    results <- caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              archsize,
              0,
              prec)
    expect_true(results$success==FALSE)
})

