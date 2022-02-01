#' Successive Pareto fronts of a population
#'
#' calculates the successive Pareto fronts of a population (classification "onion peel"), when objectives need to be maximized.
#'
#' @param matobj : matrix [ NInd , NObj ] of objectives
#' @return f : vector of dimension NInd of dominances
#'
#' @examples
#' # Definition of the parameters
#' matobj <- matrix(runif(200), 100, 2)
#' # Call the function
#' pareto_rank <- dominate(matobj)
#'
#' @author Fabrice Zaoui

dominate <- function(matobj) {
  # Test for possible NA values
  try(if(any(is.na(matobj))){
    message("The dominance cannot be calcultated due to a NA value")
    invokeRestart("abort")
  })
  # Call "pareto" function
  nind <- dim(matobj)[1] # Number of individuals
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
