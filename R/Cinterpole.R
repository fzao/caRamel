#' Cinterpole
#' 
#' La fonction Cinterpole propose n nouveaux candidats par interpolation dans des simplexes de l'espace des objectifs.
#'  
#' @param param matrice [ NPoints , NPar ] des parametres deja evalues
#' @param crit matrice [ Npoints , NObj ] des criteres associes
#' @param simplices matrice [ NSimp , NObj+1 ] contenant tout ou partie de la triangulation de l'espace des objectifs
#' @param volume matrice [ NSimp , 1 ] donnant le volume de chaque simplexe (mesure de la probabilite d'interpoler dans ce simplexe)
#' @param n nombre de nouveaux vecteurs a generer
#' @return xnouv matrice [ n , NPar ] des nouveaux vecteurs
#' @return pcrit : matrice [ n , NObj ] positions estimees des nouveaux jeux dans l'espace des objectifs
#' @author F. Zaoui
#' @export

Cinterpole <- function(param, crit, simplices, volume, n) {

  ix <- rselect(n, volume)
  nobj <- dim(crit)[2]
  nvar <- dim(param)[2]
  
  xnouv <- matrix(0, nrow = n, ncol = nvar)
  pcrit <- matrix(0, nrow = n, ncol = nobj)
  
  for (i in 1:n) {
    isimp <- ix[i]
    
    param_s <- param[simplices[isimp,],]
    crit_s <- crit[simplices[isimp,],]
    
    # Tirage aleatoire d'un jeu de coordonnees barycentriques
    xb <- runif(nobj + 1)
    xb <- xb / sum(xb)
    
    xnouv[i,] <- xb %*% param_s # Coordonnees de la C.L. dans l'espace des parameres
    pcrit[i,] <- xb %*% crit_s  # Coordonnees de la C.L. dans l'espace des criteres
  }
  return(list("xnouv" = xnouv, "pcrit" = pcrit))
}
