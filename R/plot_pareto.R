#' Plotting of a population of objectives and Pareto front
#' 
#' Plots graphs the population regarding each couple of objectives and emphasizes the Pareto front
#' 
#' @param MatObj : matrix of the objectives [NInd, nobj]
#'
#' @examples
#' # Definition of the population
#' Pop <- matrix(runif(300), 100, 3)
#' # Call the function
#' plot_pareto(MatObj = Pop)
#' 
#' @author C. Monteil

plot_pareto <- function(MatObj,nobj=NULL, objnames=NULL){
  
  if (is.null(nobj)) nobj <- length(MatObj[1,])
  if (is.null(objnames)){objnames=paste("Obj",as.character(c(1:nobj)),sep="")}
  is_pareto <- pareto(MatObj)
  
  # Graphs
  nb_fen = choose(n = nobj, k=2)
  if (nb_fen <= 4){
    par(mfrow=c(ceiling(nobj/2),2-(nb_fen==1)))
  } else {
    par(mfrow = c(3,ceiling(nb_fen/3)))
  }
  l = seq(1,nobj)
  for (i_fig in 1:(nobj-1)){
    l_tmp = l[-i_fig]; l_tmp = l_tmp[l_tmp>i_fig]
    for (i_fig2 in l_tmp){
      plot(MatObj[, i_fig], MatObj[, i_fig2],xlab = objnames[i_fig],ylab = objnames[i_fig2])
      points(MatObj[is_pareto, i_fig], MatObj[is_pareto, i_fig2],col=2,pch=16,cex=1.2)
    }
  }
  
}
