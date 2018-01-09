#' val2rank
#' 
#' La fonction val2rank convertit les valeurs d'un vecteur en leur rang 
#'  
#' @param X vecteur a traiter
#' @param opt entier qui donne la regle a suivre en cas  de rangs lies (valeurs repetees) : si opt = 1, on renvoie le rang moyen, si opt = 2, on renvoie le rang correspondant dans la serie des valeurs uniques, si opt = 3, on renvoie le rang max
#' @return R vecteur des rangs
#' @author F. Zaoui
#' @export
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
