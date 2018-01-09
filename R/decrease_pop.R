#' decrease_pop
#' 
#' La fonction decrease_pop reduit la population de jeux de parametres.
#'  
#' @param matobj matrice des objectifs, de dimension (njeux, nobj)
#' @param minmax : vecteur de booleens, de dimension nobj : TRUE si maximisation de l'objectif, FALSE sinon
#' @param prec : vecteur de dimension nobj : precision
#' @param archsize : entier, taille de l'archive
#' @param popsize : entier : taille de la population
##' @return
##' Une liste contenant deux elements:
##' \describe{
##' \item{ind_arch}{}
##' \item{ind_pop}{} 
##' }
#' @author F. Zaoui
#' @export
decrease_pop <- function(matobj, minmax, prec, archsize, popsize) {

  nobj <- dim(matobj)[2]
  ind_pop <- matrix(data = 1:dim(matobj)[1])
  ind_arch <- NULL

  matobj[, !minmax] <- -matobj[, !minmax]

  Fo <- dominate(matobj)

  # Choix des points conserves
  indices <- downsize(matobj, Fo, prec)
  pop <- matobj[indices, ]

  ind_pop <- ind_pop[indices, ]

  #////////////////////// Recalcul des fronts ////////////////////////
  Fo <- dominate(pop)

  # Tri par front croissant
  Fs <- sort(Fo, index.return = TRUE)
  pop <- pop[Fs$ix, ]
  Fo <- Fo[Fs$ix]

  ind_pop <- ind_pop[Fs$ix]

  # Separation "elite" / reste de la population
  arch <- matrix(pop[Fo == 1, ], ncol=nobj)
  pop <- matrix(pop[Fo > 1, ], ncol=nobj)
  ind_arch <- ind_pop[Fo == 1]
  ind_pop <- ind_pop[Fo > 1]

  Fo <- Fo[Fo > 1]

  #////////////////////// Reduction de la taille de la population si > popsize //////////////////////
  
  # Si la taille de la population depasse popsize, on elimine les jeux de parametres ayant un Fo eleve
  if (dim(pop)[1] > popsize) {
    Fmax <- Fo[popsize]
    ind_pop1 <- ind_pop[Fo < Fmax]

    n2 <- popsize - length(ind_pop1)

    if (n2 > 0) { # S'il reste de la place on selectionne n2 points en plus
      ix2 <- which(Fo == Fmax)
      r <- rselect(n2, matrix(1, nrow = length(ix2), ncol = 1))
      ix2 <- ix2[r]
      ind_pop1 <- c(ind_pop1, ind_pop[ix2])
    }
    ind_pop <- ind_pop1
  }

  # /////////////////// Reduction de la taille de l'archive si > archsize ///////////////
  if (dim(arch)[1] > archsize) {
    arch_down <- TRUE
    arch_prec <- 2 * prec
    while (arch_down) {
      indices <-
        downsize(arch, matrix(1, nrow = dim(arch)[1], ncol = 1), arch_prec)
      arch_down <- length(indices) > archsize
      arch_prec <- 2 * arch_prec
    }
    ind_arch <- ind_arch[indices]
  }

  return(list("arch" = ind_arch, "pop" = ind_pop))
}
