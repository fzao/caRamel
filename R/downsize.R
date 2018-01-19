#' downsize
#' 
#' La fonction downsize renvoie les indentifiants d'une population en ne conservant qu'un seul point par boite (de resolution prec)
#'  
#' @param points : matrix of objectives
#' @param Fo : rank on the front of each point (1: dominates on the Pareto)
#' @param prec : (double, length = nobj) desired accuracy for sorting objectives
#' @return vector indices
#' @author Fabrice Zaoui
#' @export

downsize <- function(points, Fo, prec) {
  the_box <- boxes(points, prec)
  
  imax1 <- apply(points, 2, which.max)
  imax1 <- sort(unique(imax1))
  
  # Conservation of a single point per box at most
  u <- sort(unique(the_box))
  indices <- NULL
  
  for (i in 1:length(u)) {
    iu <- which(the_box == u[i])
    
    # The candidates to represent the box are limited to those of smaller front
    Fmin <- min(Fo[iu])
    iux <- which(Fo[iu] == Fmin)
    iu <- iu[iux]
    iux <- NULL
    
    #If in "or" there is one (or more) set that maximizes absolutely one of the obj,
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
