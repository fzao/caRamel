#' Dimprove
#' 
#' La fonction Dimprove determine des directions d'amelioration
#'  
#' @param o_splx matrices des objectifs des simplexes (nrow = npoints, ncol = nobj)
#' @param f_splx vecteur (npoints) des numeros de Pareto associes (1 = domine)
#' @return liste des elements "oriedge" : aretes orientees et "ledge" : longueur
#' @author F. Zaoui
#' @export

Dimprove <- function(o_splx, f_splx) {
  if1 <- which(f_splx == 1)
  ns <- dim(o_splx)[1]
  nobj <- dim(o_splx)[2]
  oriedge <- NULL
  ledge <- NULL
  
  for (k in 1:length(if1)) {
    is_dominated <- dominated(o_splx[if1[k], ], o_splx)
    ideb <- which(is_dominated == TRUE)
    if (length(ideb) > 0) {
      ak <-
        cbind(matrix(ideb), matrix(if1[k], nrow = length(ideb), ncol = 1))
      dobj <- matrix(o_splx[ak[, 2], ] - o_splx[ak[, 1], ], ncol = nobj)
      lk <- matrix(sqrt(apply(dobj**2, 1, sum)), ncol = 1)
      oriedge <- rbind(oriedge, ak)
      ledge <- rbind(ledge, lk)
    }
  }
  
  return(list("oriedge" = oriedge, "ledge" = ledge))
}
