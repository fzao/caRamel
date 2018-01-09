#' dominate
#' 
#' La fonction dominate calcule les fronts de Pareto successifs d'une population (classement "en pelure d'oignon")
#'  
#' @param matobj matrice [ NInd , NObj ] des objectifs
#' @return f vecteur de dimension NInd des dominances
#' @author F. Zaoui
#' @export

dominate <- function(matobj) {
  # Appelle la fonction "pareto"
  
  nind <- dim(matobj)[1] # Nombre d'individus
  nobj <- dim(matobj)[2] # Nombre de criteres
  f <- matrix(data = 0., nrow = nind)
  ix <- which(f == 0)
  k <- 1
  while (length(ix) > 0) {
    ftmp <- pareto(matobj[ix, , drop = FALSE])
    f[ix] <- k * ftmp
    ix <- which(f == 0)
    k <- k + 1
  }
  return(f)
}
