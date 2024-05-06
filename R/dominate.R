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
#' @author Alban de Lavenne, Fabrice Zaoui

dominate <- function(matobj) {

  # Test for possible NA values
  try(if(any(is.na(matobj))){
    message("The dominance cannot be calcultated due to a NA value")
    invokeRestart("abort")
  })

  if (!is.double(matobj)) {storage.mode(matobj) <- 'double'}
  nind <- dim(matobj)[1] # Number of individuals
  ord <- order(matobj[,1], matobj[,2], decreasing=TRUE) # order to make it faster
  res <- matrix(data = .Call(c_dominate, matobj[ord,]), nrow = nind)
  res[ord] <- res # original order
  return(res)

}
