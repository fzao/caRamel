#' Cextrap
#' 
#' gives n new candidates by extrapolation along orthogonal directions to the Pareto front in the space of the objectives
#'  
#' @param param : matrix [ NPoints , NPar ] of already evaluated parameters
#' @param crit : matrix [ Npoints , NObj ] of associated criteria
#' @param directions : matrix [ NDir, 2 ] the starting and ending points of the candidate vectors
#' @param longu : matrix [ NDir , 1 ] giving the length of each segment thus defined in the OBJ space (measure of the probability of exploring this direction)
#' @param n : number of new vectors to generate
#' @return xnouv : matrix [ n , NPar ] of new vectors
#' @return pcrit : matrix [ n , NObj ] estimated positions of new sets in the goal space
#' @author Fabrice Zaoui

Cextrap <- function(param, crit, directions, longu, n) {
  
  ix <- rselect(n, longu)
  nobj <- dim(crit)[2]
  nvar <- dim(param)[2]
  
  xnouv <- matrix(0, nrow = n, ncol = nvar)
  pcrit <- matrix(0, nrow = n, ncol = nobj)
  
  for (i in 1:n) {
    iar <- ix[i]
    
    param_a <- as.matrix(param[directions[iar, ], ])
    crit_a <- crit[directions[iar, ], ]
    
    # Search vector
    u_param <- diff(param_a)
    u_crit <- diff(crit_a)
    
    boost <- mean(longu) / longu[iar] # The shorter the edge, the more "boost" the search
    lambda <- -log(1 - runif(1))
    
    xnouv[i, ] <- param_a[2, ] + boost * lambda * u_param # Coordinates in the space of the parameters
    pcrit[i, ] <- crit_a[2, ]  + boost * lambda * u_crit  # Coordinates in the space of the criteria
  }
  return(list("xnouv" = xnouv, "pcrit" = pcrit))
}
