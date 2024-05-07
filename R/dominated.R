#'Rows domination of a matrix by a vector
#'
#' indicates which rows of the matrix Y are dominated by the vector (row) x
#'
#' @param x : row vecteur
#' @param Y : matrix
#' @return D : vector of booleans
#'
#' @examples
#' # Definition of the parameters
#' Y <- matrix(rexp(200), 100, 2)
#' x <- Y[1,]
#' # Call the function
#' res <- dominated(x, Y)
#'
#' @author Alban de Lavenne, Fabrice Zaoui

dominated <- function(x, Y) {

  if (!is.double(x)) {storage.mode(x) <- 'double'}
  if (!is.double(Y)) {storage.mode(Y) <- 'double'}
  res <- .Call(c_dominated, x, Y)
  return(res == 1)

}
