#' downsize
#' 
#' La fonction downsize renvoie les indentifiants d'une population en ne conservant qu'un seul point par boite (de resolution prec)
#'  
#' @param points : matrice des objectifs
#' @param Fo : rang sur le front de chaque point (1 : domine sur le Pareto)
#' @param prec (double, length = nobj) precision souhaitee pour le tri des objectifs
#' @return vecteur des indices
#' @author F. Zaoui
#' @export

downsize <- function(points, Fo, prec) {
  the_box <- boxes(points, prec)
  
  imax1 <- apply(points, 2, which.max)
  imax1 <- sort(unique(imax1))
  
  # Conservation d'un seul point par boite au maximum
  u <- sort(unique(the_box))
  indices <- NULL
  
  for (i in 1:length(u)) {
    iu <- which(the_box == u[i])
    
    # Les candidats a representer la boite sont limites a ceux de plus petit front
    Fmin <- min(Fo[iu])
    iux <- which(Fo[iu] == Fmin)
    iu <- iu[iux]
    iux <- NULL
    
    #Si dans "iu" il y a un (ou plusieurs) jeu qui maximise(nt) absolument l'un des obj, c'est lui qu'on garde
    # (cas ou la boite contient au moins un point du front)
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
        # Sinon, Tirage aleatoire d'un nombre entre 1 et length(iu)
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
