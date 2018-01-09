#' Crecombination
#' 
#' La fonction Crecombination effectue une recombinaison des jeux de parametres.
#'  
#' @param param matrice [ . , NPar ] de la population des parametre
#' @param blocs liste de vecteurs d'entier : liste de blocs de variables pour la recombinaison
#' @param n nombre de nouveaux vecteurs a generer
#' @return xnew matrice [ n , NPar ] des nouveaux vecteurs
#' @author F. Zaoui
#' @export

Crecombination <- function(param, blocs, n) {
  npar <- dim(param)[2]
  paramiso <- matrix(TRUE, ncol = npar)
  nblocs <- length(blocs)
  if (nblocs > 0) {
    for (i in 1:nblocs) {
      paramiso[blocs[[i]]] <- FALSE
    }
  }
  ix <- which(paramiso)
  
  blocs1 <- blocs
  n1 <- nblocs + length(ix)
  
  for (i in 1:length(ix)) {
    blocs1[[(nblocs+i)]] <- c(ix[i])
  }
  
  xnouv <- matrix(0, nrow = n, ncol = npar)
  for (i in 1:n) {
    gamebloc <- rselect(n1, matrix(1, nrow = dim(param)[1], ncol = 1)) # On selectionne n1 jeux dont on prendra n1 blocs
    p <- matrix(0, ncol = npar, nrow = 1)
    for (b in 1:n1) {
      p[blocs1[[b]]] <- param[gamebloc[b], blocs1[[b]]]
    }
    xnouv[i,] <- p
  }

  return(xnouv)
}
