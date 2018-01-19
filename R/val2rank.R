#' val2rank
#' 
#' converts the values of a vector into their rank
#'  
#' @param X : vector to treat
#' @param opt : integer which gives the rule to follow in case of ranks linked (repeated values): if opt = 1, one returns the average rank, if opt = 2, one returns the corresponding rank in the series of the unique values, if opt = 3, return the max rank
#' @return R : rank vector
#' @author Fabrice Zaoui

val2rank <- function(X, opt) {

  Rs <- matrix(0., nrow = length(X), ncol = 1)
  R <- matrix(0., nrow = length(X), ncol = 1)
  
  xso <- sort(X, index.return = TRUE)
  Xs <- xso$x
  ix <- xso$ix
  d <- c(diff(Xs), 1)
  ideb <- 1
  u <- 0
  
  for (i in 1:length(Xs)) {
    if (d[i] != 0) {
      if (opt == 1) {
        Rs[ideb:i, 1] <- (ideb + i) / 2
      } else if (opt == 2) {
        u <- u + 1
        Rs[ideb:i, 1] <- u
      } else{
        Rs[ideb:i, 1] <- i
      }
      ideb <- i + 1
    }
  }
  R[ix, ] <- Rs
  return(R)
}
