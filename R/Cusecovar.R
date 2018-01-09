#' Cusecovar
#' 
#' La fonction Cusecovar Propose de nouveaux vecteurs de parametres respectant une structure de covariance.
#'  
#' @param xref matrice [ . , NPar ] de la population de reference dont on souhaite utiliser la structure de covariance
#' @param amplif facteur d'amplification de l'ecart-type sur chaque parametre
#' @param n nombre de nouveaux vecteurs a generer
#' @return xnew matrice [ n , NPar ] des nouveaux vecteurs
#' @author F. Zaoui
#' @export

Cusecovar <- function(xref, amplif, n) {
  
  g <- apply(xref, 2, mean) # barycentre de la population de reference (dans l'espace des parametres)
  
  # Calcul de la matrice de variances-covariances des lignes
  rr <- matvcov(xref, g)
  
  # Selection des parametres de variance non-nulle et amplification
  sref <- apply(xref, 2, sd)
  sref[which(sref <= 0)] <- NaN
  pvar <- which(!is.na(sref))
  rr <- rr[pvar, pvar]
  d <- amplif * diag(sref[pvar], length(pvar))

  # Calcul de la matrice de Cholesky qui va bien
  tc <- chol(d %*% rr %*% d)
  
  # Initialisation des nouveaux jeux au barycentre
  xnew <- matrix(g,
                 nrow = n,
                 ncol = length(g),
                 byrow = TRUE)
  
  # Tirage aleatoire du bruit
  eps <- matrix(rnorm(n * length(pvar)), nrow = n)
  
  # Superposition
  xnew[, pvar] <- xnew[, pvar] + eps %*% tc
  
  return(xnew)
}
