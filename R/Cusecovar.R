#' Cusecovar
#' 
#' proposes new parameter vectors respecting a covariance structure
#'  
#' @param xref : matrix [ . , NPar ] of the reference population whose covariance structure is to be used
#' @param amplif : amplification factor of the standard deviation on each parameter
#' @param n : number of new vectors to generate
#' @return xnew : matrix [ n , NPar ] of new vectors
#' 
#' @examples
#' # Definition of the parameters
#' xref <- matrix(rexp(35), 35, 1)
#' amplif <- 2.
#' n <- 5
#' # Call the function
#' res <- Cusecovar(xref, amplif, n)
#' 
#' @author Fabrice Zaoui

Cusecovar <- function(xref, amplif, n) {
  
  g <- apply(xref, 2, mean) # center of reference population (in the parameter space)

  # Calculation of the variances-covariances matrix of the lines
  rr <- matvcov(xref, g)
  
  # Selection of parameters of non-zero variance and amplification
  sref <- apply(xref, 2, sd)
  sref[which(sref <= 0)] <- NaN
  pvar <- which(!is.na(sref))
  rr <- rr[pvar, pvar]
  d <- amplif * diag(sref[pvar], length(pvar))

  # Calculation of the Cholesky matrix
  tc <- chol(d %*% rr %*% d)
  
  # Initialization of new sets at the center of gravity
  xnew <- matrix(g,
                 nrow = n,
                 ncol = length(g),
                 byrow = TRUE)
  
  # Random draw of noise
  eps <- matrix(rnorm(n * length(pvar)), nrow = n)
  
  # Superposition
  xnew[, pvar] <- xnew[, pvar] + eps %*% tc
  
  return(xnew)
}
