#' Downsizing of a population to only one individual per box up to a given accuracy
#'
#' reduces the number of individuals in a population to only one individual per box up to a given accuracy
#'
#' @param points : matrix of objectives
#' @param Fo : rank on the front of each point (1: dominates on the Pareto)
#' @param prec : (double, length = nobj) desired accuracy for sorting objectives
#' @return vector indices
#'
#' @examples
#' # Definition of the parameters
#' points <- matrix(rexp(200), 100, 2)
#' prec <- c(1.e-3, 1.e-3)
#' Fo <- sample(1:100, 100)
#' # Call the function
#' res <- downsize(points, Fo, prec)
#'
#' @author Fabrice Zaoui

downsize <- function(points, Fo, prec) {
  the_box <- boxes(points, prec)

  imax1 <- apply(points, 2, which.max)
  imax1 <- sort(unique(imax1))

  # Conservation of a single point per box at most
  u <- sort(unique(the_box))
  indices <- NULL

  for (i in seq_len(length(u))) {
    iu <- which(the_box == u[i])

    # The candidates to represent the box are limited to those of smaller front
    Fmin <- min(Fo[iu])
    iux <- which(Fo[iu] == Fmin)
    iu <- iu[iux]
    iux <- NULL

    # If in "or" there is one (or more) set that maximizes absolutely one of the obj,
    # it is him that one keeps (case where the box contains at least one point of the front)
    if (Fmin == 1) {
      D <-
        matrix(iu, nrow = length(iu), ncol = 1) %*% matrix(1, nrow = 1, ncol =
                                                             length(imax1)) -
        matrix(1, nrow = length(iu), ncol = 1) %*% matrix(imax1, nrow =
                                                            1, ncol = length(imax1))
      iux <- isTRUE(which(D == 0) > 0)
    }

    if (!is.null(iux)) {
      if (iux == TRUE) {
        iu <- iu[1]
      } else {
        # Otherwise, random draw of a number between 1 and length(iu)
        j <- ceiling(length(iu) * runif(1))
        iu <- iu[j]
      }
    } else{
      j <- ceiling(length(iu) * runif(1))
      iu <- iu[j]
    }

    indices <- c(indices, iu)
  }
  return(indices)
}
