#' Decreasing of the population of parameters sets
#'
#' decreases the population of parameters sets
#'
#' @param matobj : matrix of objectives, dimension (ngames, nobj)
#' @param minmax : vector of booleans, of dimension nobj: TRUE if maximization of the objective, FALSE otherwise
#' @param prec : nobj dimension vector: accuracy
#' @param archsize : integer: archive size
#' @param popsize : integer: population size
##' @return
##' A list containing two elements:
##' \describe{
##' \item{ind_arch}{indices of individuals in the updated Pareto front}
##' \item{ind_pop}{indices of individuals in the updated population}
##' }
#'
#' @examples
#' # Definition of the parameters
#' matobj <- matrix(rexp(200), 100, 2)
#' prec <- c(1.e-3, 1.e-3)
#' archsize <- 100
#' minmax <- c(FALSE, FALSE)
#' popsize <- 100
#' # Call the function
#' res <- decrease_pop(matobj, minmax, prec, archsize, popsize)
#'
#' @author Fabrice Zaoui

decrease_pop <- function(matobj, minmax, prec, archsize, popsize) {

  nobj <- dim(matobj)[2]
  ind_pop <- matrix(data = 1:dim(matobj)[1])
  ind_arch <- NULL

  matobj[, !minmax] <- -matobj[, !minmax]

  Fo <- dominate(matobj)

  # Choice of retained points
  indices <- downsize(matobj, Fo, prec)
  pop <- matobj[indices, ]
  dim(pop) <- c(length(indices), nobj)

  ind_pop <- ind_pop[indices, ]

  #////////////////////// Recalculation of fronts ////////////////////////
  Fo <- dominate(pop)

  # Sort by increasing front
  Fs <- sort(Fo, index.return = TRUE)
  pop <- pop[Fs$ix, ]
  dim(pop) <- c(length(Fs$ix), nobj)
  Fo <- Fo[Fs$ix]

  ind_pop <- ind_pop[Fs$ix]

  # Separation "elite" / rest of the population
  arch <- matrix(pop[Fo == 1, ], ncol = nobj)
  pop <- matrix(pop[Fo > 1, ], ncol = nobj)
  ind_arch <- ind_pop[Fo == 1]
  ind_pop <- ind_pop[Fo > 1]

  Fo <- Fo[Fo > 1]

  # Decreasing the population size
  # If the population size exceeds popsize, elimination of the parameter sets with a high Fo
  if (dim(pop)[1] > popsize) {
    Fmax <- Fo[popsize]
    ind_pop1 <- ind_pop[Fo < Fmax]

    n2 <- popsize - length(ind_pop1)

    if (n2 > 0) { # If still possible, selection of n2 points in addition
      ix2 <- which(Fo == Fmax)
      r <- rselect(n2, matrix(1, nrow = length(ix2), ncol = 1))
      ix2 <- ix2[r]
      ind_pop1 <- c(ind_pop1, ind_pop[ix2])
    }
    ind_pop <- ind_pop1
  }

  #  Decreasing the archive size
  if (dim(arch)[1] > archsize) {
    arch_down <- TRUE
    arch_prec <- 2 * prec
    while (arch_down) {
      indices <-
        downsize(arch, matrix(1, nrow = dim(arch)[1], ncol = 1), arch_prec)
      arch_down <- length(indices) > archsize
      arch_prec <- 2 * arch_prec
    }
    ind_arch <- ind_arch[indices]
  }

  return(list("arch" = ind_arch, "pop" = ind_pop))
}
