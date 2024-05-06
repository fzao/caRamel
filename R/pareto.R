#' Indicates which rows are Pareto
#'
#' indicates which rows of the X criterion matrix are Pareto, when objectives need to be maximized
#'
#' @param X : matrix of objectives [NInd * NObj]
#' @return Ft : vector [NInd], TRUE when the set is on the Pareto front.
#'
#' @examples
#' # Definition of the parameters
#' X <- matrix(runif(200), 100, 2)
#' # Call the function
#' is_pareto <- pareto(X)
#'
#' @author Alban de Lavenne, Fabrice Zaoui

pareto <- function(X) {

  if (!is.double(X)) {storage.mode(X) <- 'double'}
  ord <- order(X[,1],X[,2],decreasing=TRUE) # order to make it faster
  res <- .Call(c_pareto, X[ord,]) > 0
  res[ord] <- res # original order
  return(res)

}
