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
#' @author Fabrice Zaoui

pareto <- function(X) {
  
  Ft <- logical(length = dim(X)[1]) # Initially all vectors are marked as non-Pareto
  
  Xtmp <- X
  itmp <- matrix(data = 1:dim(X)[1])
  ix1 <- 1
  i1 <- 1
  
  while (1) {
    X1 <- Xtmp[ix1, ]
    is_dominated <- dominated(X1, Xtmp)
    Xtmp <- Xtmp[!is_dominated, ]
    itmp <- itmp[!is_dominated]
    ix1 <- which(itmp > i1)[1]
    if (is.na(ix1))
      break
    i1 <- itmp[ix1]
  }
  
  Ft[itmp] <- TRUE
  return(Ft)
  
}
