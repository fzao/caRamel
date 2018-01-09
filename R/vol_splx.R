#' vol_splx
#' 
#' La fonction vol_splx calcule le volume d'un simplex
#'  
#' @param S : matrice (d+1) lignes * d colonnes contenant les coordonnees en d-dim des d+1 sommets d'un simplex
#' @return V volume du simplex
#' @author F. Zaoui
#' @export

vol_splx <- function(S) {
  
  d <- dim(S)[2]
  A <- S[2:dim(S)[1], ] - matrix(1, nrow = d, ncol = 1) %*% S[1, ]
  V <- (1 / factorial(d)) * abs(det(A))
  
  return(V)
}
