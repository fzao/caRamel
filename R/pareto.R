#' pareto
#' 
#' La fonction pareto indique quelles lignes de la matrice de criteres X sont Pareto
#'  
#' @param X matrice [NInd * NObj]
#' @return Ft matrice colonne [NInd * 1]
#' @author F. Zaoui
#' @export

pareto <- function(X) {
  
  Ft <- logical(length = dim(X)[1]) # Initialement tous les vecteurs sont marques comme non-Pareto
  
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
