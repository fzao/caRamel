#' pareto
#' 
#' indicates which rows of the X criterion matrix are Pareto
#'  
#' @param X : matrix [NInd * NObj]
#' @return Ft : column matrix [NInd * 1]
#' 
#' @examples
#' # Definition of the parameters
#' X <- matrix(rexp(200), 100, 2)
#' # Call the function
#' res <- dominate(X)
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
