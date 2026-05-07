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
  if (!is.double(X)) storage.mode(X) <- 'double'
  nobj <- ncol(X)
  
  if (nobj == 2L) {
    ord <- order(X[, 1], X[, 2], decreasing = TRUE)
    res <- .Call(c_pareto_2d, X[ord, , drop = FALSE]) > 0
  } else if (nobj == 3L) {
    ord <- order(X[, 1], X[, 2], X[, 3], decreasing = TRUE)
    res <- .Call(c_pareto_3d, X[ord, , drop = FALSE]) > 0
  } else {
    ord <- do.call(order, c(lapply(seq_len(nobj), function(k) X[, k]),
                            list(decreasing = TRUE)))
    res <- .Call(c_pareto, X[ord, , drop = FALSE]) > 0
  }
  out <- logical(length(res))
  out[ord] <- res
  out
}
