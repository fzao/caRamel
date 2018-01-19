#' Cinterp
#' 
#' proposes n new candidates by interpolation in simplexes of the objective space
#'  
#' @param param : matrix [ NPoints , NPar ] of already evaluated parameters
#' @param crit : matrix [ Npoints , NObj ] of associated criteria
#' @param simplices : matrix [ NSimp , NObj+1 ] containing all or part of the triangulation of the space of the objectives
#' @param volume : matrix [ NSimp , 1 ] giving the volume of each simplex (measure of the probability of interpolating in this simplex)
#' @param n : number of new vectors to generate
#' @return xnouv : matrix [ n , NPar ] of new vectors
#' @return pcrit : matrix [ n , NObj ] estimated positions of new sets in the goal space
#' @author Fabrice Zaoui
#' @export

Cinterp <- function(param, crit, simplices, volume, n) {

  ix <- rselect(n, volume)
  nobj <- dim(crit)[2]
  nvar <- dim(param)[2]
  
  xnouv <- matrix(0, nrow = n, ncol = nvar)
  pcrit <- matrix(0, nrow = n, ncol = nobj)
  
  for (i in 1:n) {
    isimp <- ix[i]
    
    param_s <- param[simplices[isimp,],]
    crit_s <- crit[simplices[isimp,],]
    
    # Random drawing of a set of barycentric coordinates
    xb <- runif(nobj + 1)
    xb <- xb / sum(xb)
    
    xnouv[i,] <- xb %*% param_s # Coordinates in the space of the parameters
    pcrit[i,] <- xb %*% crit_s  # Coordinates in the space of the criteria
  }
  return(list("xnouv" = xnouv, "pcrit" = pcrit))
}
