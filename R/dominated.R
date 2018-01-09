#' dominated
#' 
#' La fonction dominated indique quelles lignes de la matrice Y sont dominees par le vecteur (ligne) x
#'  
#' @param x : vecteur ligne
#' @param Y : matrice
#' @return D : vecteur de booleens
#' @author F. Zaoui
#' @export

dominated <- function(x, Y) {
  
  # Fonction predicat indiquant quelles lignes de la matrice Y
  # sont dominees par le vecteur (ligne) x
  
  X = matrix(rep(x, dim(Y)[1]), ncol = dim(Y)[2], byrow = TRUE)
  
  D = X - Y
  nobj <- dim(D)[2]
  res1 <- rowSums(D >= 0) == nobj
  res2 <- rowSums(D > 0) > 0
  
  return(res1 & res2)
}
