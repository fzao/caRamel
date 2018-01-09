#' Cextrapole
#' 
#' La fonction Cextrapole propose n nouveaux candidats par extrapolation le long de directions orthogonales au Front de Pareto dans l'espace des objectifs
#'  
#' @param param matrice [ NPoints , NPar ] des parametres deja evalues
#' @param crit matrice [ Npoints , NObj ] des criteres associes
#' @param directions  matrice [ NDir, 2 ] les points de debut et de fin des vecteurs directeurs candidats
#' @param longu matrice [ NDir , 1 ] donnant la longueur de chaque segment ainsi defini dans l'espace des OBJ (mesure de la probabilite d'explorer cette direction)
#' @param n nombre de nouveaux vecteurs a generer
#' @return xnouv matrice [ n , NPar ] des nouveaux vecteurs
#' @return pcrit : matrice [ n , NObj ] positions estimees des nouveaux jeux dans l'espace des objectifs
#' @author F. Zaoui
#' @export


Cextrapole <- function(param, crit, directions, longu, n) {
  
  ix <- rselect(n, longu)
  nobj <- dim(crit)[2]
  nvar <- dim(param)[2]
  
  xnouv <- matrix(0, nrow = n, ncol = nvar)
  pcrit <- matrix(0, nrow = n, ncol = nobj)
  
  for (i in 1:n) {
    iar <- ix[i]
    
    param_a <- as.matrix(param[directions[iar, ], ])
    crit_a <- crit[directions[iar, ], ]
    
    # Vecteur de recherche
    u_param <- diff(param_a)
    u_crit <- diff(crit_a)
    
    boost <- mean(longu) / longu[iar] # Plus l'arete est courte plus on "booste" la recherche
    lambda <- -log(1 - runif(1))
    
    xnouv[i, ] <- param_a[2, ] + boost * lambda * u_param # Coordonnees de la C.L. dans l'espace des parametres
    pcrit[i, ] <- crit_a[2, ]  + boost * lambda * u_crit  # Coordonnees de la C.L. dans l'espace des criteres
  }
  return(list("xnouv" = xnouv, "pcrit" = pcrit))
}
