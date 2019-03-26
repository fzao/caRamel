#' Recombination of the sets of parameters
#' 
#' performs a recombination of the sets of parameters
#'  
#' @param param : matrix [ . , NPar ] of the population of parameters
#' @param blocks : list of integer vectors: list of variable blocks for recombination
#' @param n : number of new vectors to generate
#' @return xnew : matrix [ n , NPar ] of new vectors
#' 
#' @examples
#' # Definition of the parameters
#' param <- matrix(rexp(15), 15, 1)
#' blocks <- NULL
#' n <- 5
#' # Call the function
#' res <- Crecombination(param, blocks, n)
#' 
#' @author Fabrice Zaoui

Crecombination <- function(param, blocks, n) {

  npar <- dim(param)[2]
  paramiso <- matrix(TRUE, ncol = npar)
  nblocks <- length(blocks)
  if (nblocks > 0) {
    for (i in 1:nblocks) {
      paramiso[blocks[[i]]] <- FALSE
    }
  }
  ix <- which(paramiso)
  
  blocks1 <- blocks
  n1 <- nblocks + length(ix)
  
  for (i in 1:length(ix)) {
    blocks1[[(nblocks+i)]] <- c(ix[i])
  }
  
  xnew <- matrix(0, nrow = n, ncol = npar)
  for (i in 1:n) {
    gamebloc <- rselect(n1, matrix(1, nrow = dim(param)[1], ncol = 1)) # Selection of n1 sets for n1 blocks
    p <- matrix(0, ncol = npar, nrow = 1)
    for (b in 1:n1) {
      p[blocks1[[b]]] <- param[gamebloc[b], blocks1[[b]]]
    }
    xnew[i,] <- p
  }

  return(xnew)
}
