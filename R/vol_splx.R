#' vol_splx
#' 
#' calculates the volume of a simplex
#'  
#' @param S : matrix (d+1) rows * d columns containing the coordinates in d-dim of d + 1 vertices of a simplex
#' @return V : simplex volume
#' @author Fabrice Zaoui
#' @export

vol_splx <- function(S) {
  
  d <- dim(S)[2]
  A <- S[2:dim(S)[1], ] - matrix(1, nrow = d, ncol = 1) %*% S[1, ]
  V <- (1 / factorial(d)) * abs(det(A))
  
  return(V)
}
