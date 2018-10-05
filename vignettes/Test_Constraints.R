## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----caRa----------------------------------------------------------------
library(caRamel)

## ----constr--------------------------------------------------------------
constr_ex <- function(i) {
  # functions f1 and f2
  s1 <- x[i,1]
  s2 <- (1. + x[i,2]) / x[i,1]
  # now test for the feasibility
  # constraint g1
  if((x[i,2] + 9. * x[i,1] - 6.) < 0. | (-x[i,2] + 9. * x[i,1] -1.) < 0.) {
    s1 <- NaN
    s2 <- NaN
  }
  return(c(s1, s2))
}

## ----constr_variable-----------------------------------------------------
nvar <- 2 # number of variables
bounds <- matrix(data = 0., nrow = nvar, ncol = 2) # upper and lower bounds
bounds[1, 1] <- 0.1
bounds[1, 2] <- 1.
bounds[2, 1] <- 0.
bounds[2, 2] <- 5.

## ----constr_objectives---------------------------------------------------
nobj <- 2 # number of objectives
minmax <- c(FALSE, FALSE) # min and min

## ----constr_param--------------------------------------------------------
popsize <- 100 # size of the genetic population
archsize <- 100 # size of the archive for the Pareto front
maxrun <- 1000 # maximum number of calls
prec <- matrix(1.e-3, nrow = 1, ncol = nobj) # accuracy for the convergence phase

## ----schaffer_launch, fig.show="hide", results="hide"--------------------
results <-
  caRamel(nobj,
          nvar,
          minmax,
          bounds,
          constr_ex,
          popsize,
          archsize,
          maxrun,
          prec,
          carallel=FALSE) # no parallelism

## ----schaffer_OK---------------------------------------------------------
print(results$success==TRUE)

## ----schaffer_plot1------------------------------------------------------
plot(results$objectives[,1], results$objectives[,2], main="Constr_Ex Pareto front", xlab="Objective #1", ylab="Objective #2")

## ----schaffer_plot2------------------------------------------------------
plot(results$parameters, main="Corresponding values for X", xlab="Element of the archive", ylab="X Variable")

