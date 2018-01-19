#' dominated
#' 
#' Tindicates which rows of the matrix Y are dominated by the vector (row) x
#'  
#' @param x : row vecteur
#' @param Y : matrix
#' @return D : vector of booleans
#' @author F. Zaoui

dominated <- function(x, Y) {
  
  # Predicate function indicating which rows of the matrix Y are dominated by the vector (row) x
  X = matrix(rep(x, dim(Y)[1]), ncol = dim(Y)[2], byrow = TRUE)
  
  D = X - Y
  nobj <- dim(D)[2]
  res1 <- rowSums(D >= 0) == nobj
  res2 <- rowSums(D > 0) > 0
  
  return(res1 & res2)
}
