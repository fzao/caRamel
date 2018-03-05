#' dominate
#' 
#' calculates the successive Pareto fronts of a population (classification "onion peel")
#'  
#' @param matobj : matrix [ NInd , NObj ] of objectives
#' @return f : vector of dimension NInd of dominances
#' 
#' @examples
#' # Definition of the parameters
#' matobj <- matrix(rexp(200), 100, 2)
#' # Call the function
#' res <- dominate(matobj)
#' 
#' @author Fabrice Zaoui

dominate <- function(matobj) {
  # Call "pareto" function
  nind <- dim(matobj)[1] # Number of individuals
  nobj <- dim(matobj)[2] # Number of criteria
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
