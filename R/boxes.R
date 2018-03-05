#' boxes
#' 
#' This function returns a box number for each points individual of the population
#'  
#' @param points : matrix of the objectives
#' @param prec : (double, length = nobj) desired accuracy for the objectives (edges of the boxes)
#' @return vector of numbers for the boxes. boxes[i] gives the number of the box containing points[i].
#' 
#' @examples
#' # Definition of the parameters
#' points <- matrix(rexp(200), 100, 2)
#' prec <- c(1.e-3, 1.e-3)
#' # Call the function
#' res <- boxes(points, prec)
#' 
#' @author Fabrice Zaoui

boxes <- function(points, prec) {
  
  n <- dim(points)[1]
  d <- dim(points)[2]

  valmin <- apply(points,2,min)
  val <- matrix(rep(valmin,n), ncol=d, byrow=TRUE)
  
  the_box <- points - val
  the_prec <- matrix(rep(prec,n), ncol=d, byrow=TRUE)
  the_box <- floor(the_box / the_prec)
  
  for(i in 1:d) {
    C <- val2rank(the_box[,i],2)
    the_box[,i] <- C
  }
  
  base <- 1 + apply(the_box,2,max)
  base <- c(1, base[1:length(base)-1])
  base <- cumprod(base)
  base <- matrix(rep(base,n),ncol=d, byrow=TRUE)

  return(rowSums(the_box*base))
}