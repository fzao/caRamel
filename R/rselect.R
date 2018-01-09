#' rselect
#' 
#' La fonction rselect effectue une selection de n points dans facp
#'  
#' @param n : nombre de points a selectionner
#' @param facp : vecteur de points initiaux
#' @return ix : rangs des points selectionnes (vecteur de dimension n)
#' @author F. Zaoui
#' @export

rselect <- function(n, facp) {
  if (sum(facp) > 0 & n > 0) {
    facp1 <- facp / sum(facp)
    Fc <- cumsum(facp1)
    
    ix <- matrix(0, nrow = n, ncol = 1)
    
    for (i in 1:n) {
      f <- runif(1)
      b <- which(Fc > f)
      ix[i] <- b[1]
    }
  } else{
    ix <- NULL
  }
  return(ix)
}
