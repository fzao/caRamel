#
# Test the optimization process
#

context("success")

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
archsize <- 10 # size of the archive for the Pareto front
maxrun <- 100 # maximum number of calls
prec <- matrix(1.e-3, nrow = 1, ncol = nobj) # accuracy for the convergence phase

test_that("Optimization process is OK", {
    # flag must be TRUE
    results <-
      caRamel(nobj,
              nvar,
              minmax,
              bounds,
              fonseca,
              popsize,
              archsize,
              maxrun,
              prec,
              carallel=FALSE)
    expect_true(results$success==TRUE)
})

