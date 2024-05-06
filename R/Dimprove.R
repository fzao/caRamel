#' Determination of directions for improvement
#'
#' determines directions for improvement
#'
#' @param o_splx : matrix of objectives of simplexes (nrow = npoints, ncol = nobj)
#' @param f_splx : vector (npoints) of associated Pareto numbers (1 = dominated)
#' @return list of elements "oriedge": oriented edges and "ledge": length
#'
#' @examples
#' # Definition of the parameters
#' o_splx <- matrix(rexp(6), 3, 2)
#' f_splx <- c(1,1,1)
#' # Call the function
#' res <- Dimprove(o_splx, f_splx)
#'
#' @author Fabrice Zaoui

Dimprove <- function(o_splx, f_splx) {
  if1 <- which(f_splx == 1)
  ns <- dim(o_splx)[1]
  nobj <- dim(o_splx)[2]
  oriedge <- NULL
  ledge <- NULL

  for (k in seq_len(length(if1))) {
    is_dominated <- dominated(o_splx[if1[k], ], o_splx)
    ideb <- which(is_dominated == TRUE)
    if (length(ideb) > 0) {
      ak <-
        cbind(matrix(ideb), matrix(if1[k], nrow = length(ideb), ncol = 1))
      dobj <- matrix(o_splx[ak[, 2], ] - o_splx[ak[, 1], ], ncol = nobj)
      lk <- matrix(sqrt(apply(dobj**2, 1, sum)), ncol = 1)
      oriedge <- rbind(oriedge, ak)
      ledge <- rbind(ledge, lk)
    }
  }

  return(list("oriedge" = oriedge, "ledge" = ledge))
}
