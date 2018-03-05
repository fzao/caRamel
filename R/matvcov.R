#' matvcov
#' 
#' calculates the variances-covariances matrix on the reference population
#'  
#' @param x : population
#' @param g : center of reference population (in the parameter space)
#' @return rr : variances-covariances matrix on the reference population
#' 
#' @examples
#' # Definition of the parameters
#' x <- matrix(rexp(30), 30, 1)
#' g <- mean(x)
#' # Call the function
#' res <- matvcov(x, g)
#' 
#' @author Fabrice Zaoui

matvcov <- function(x, g) {
  
  n <- dim(x)[2]
  rr <- matrix(0, nrow = n, ncol = n)

  xred <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2])
  ss <- apply(x, 2, sd)
  ss[which(ss <= 0)] <- NaN
  
  for (i in 1:n) {
    xred[, i] <- (x[, i] - g[i]) / ss[i]
  }
  
  if (n > 1) {
    for (i1 in 1:(n - 1)) {
      for (i2 in (i1 + 1):n) {
        rr[i1, i2] <- mean(xred[, i1] * xred[, i2])
      }
    }
  }
  
  rr <- rr + t(rr) + diag(1, n)
  
  return(rr)
}
