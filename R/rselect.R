#' rselect
#' 
#' performs a selection of n points in facp
#'  
#' @param n : number of points to select
#' @param facp : vector of initial points
#' @return ix : ranks of selected points (vector of dimension n)
#' @author Fabrice Zaoui
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
