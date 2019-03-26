#' Plotting of Pareto Front
#' 
#' Plot graphs of the Pareto front and a graph of optimization evolution
#' 
#' @param nobj : number of objectives
#' @param ngen : number of generations
#' @param nrun : number of model evaluations
#' @param objnames : vector of objectives names
#' @param MatObj : matrix of the objectives of the Pareto front
#' @param MatEvol : matrix of the evolution of the optimal objectives (optional)
#'
#' @examples
#' # Definition of the test function
#' viennet <- function(i) {
#'   val1 <- 0.5*(x[i,1]*x[i,1]+x[i,2]*x[i,2])+sin(x[i,1]*x[i,1]+x[i,2]*x[i,2])
#'   val2 <- 15+(x[i,1]-x[i,2]+1)*(x[i,1]-x[i,2]+1)/27+(3*x[i,1]-2*x[i,2]+4)*(3*x[i,1]-2*x[i,2]+4)/8
#'   val3 <- 1/(x[i,1]*x[i,1]+x[i,2]*x[i,2]+1) -1.1*exp(-(x[i,1]*x[i,1]+x[i,2]*x[i,2]))
#'   return(c(val1,val2,val3))
#' }
#' nobj <- 3 # Number of objectives
#' nvar <- 2 # Number of variables
#' minmax <- c(FALSE, FALSE, FALSE) # All the objectives are to be minimized
#' bounds <- matrix(data = 1, nrow = nvar, ncol = 2) # Define the bound constraints
#' bounds[, 1] <- -3 * bounds[, 1]
#' bounds[, 2] <- 3 * bounds[, 2]
#' 
#' # Caramel optimization
#' results <- caRamel(nobj, nvar, minmax, bounds, viennet, popsize = 100, archsize = 100,
#'           maxrun = 500, prec = matrix(1.e-3, nrow = 1, ncol = nobj), carallel = FALSE)
#' 
#' # Plot of results
#' ngen <- length(results$save_crit[,1])
#' nrun <- results$save_crit[ngen,1]
#' objnames <- c("Obj1","Obj2","Obj3")
#' plot_pareto(nobj, ngen, nrun, objnames, MatObj = results$objectives, MatEvol = results$save_crit)
#' 
#' @author C. Monteil

plot_pareto <- function(nobj, ngen, nrun, objnames, MatObj, MatEvol=NULL){
  
  # Graphs
  info<-paste("ngen=",ngen,", nrun=",nrun,sep="")
  nb_fen = choose(n = nobj, k=2)+1
  if (nb_fen <= 4){
    par(mfrow=c((nobj-1),2))
  } else {
    par(mfrow = c(3,floor(nb_fen/3)+1))
  }
  l = seq(1,nobj)
  for (i_fig in 1:(nobj-1)){
    l_tmp = l[-i_fig]; l_tmp = l_tmp[l_tmp>i_fig]
    for (i_fig2 in l_tmp){
      plot(MatObj[, i_fig], MatObj[, i_fig2],xlab = objnames[i_fig],ylab = objnames[i_fig2])
    }
  }
  
  if(!is.null(MatEvol)){
    plot(x=MatEvol[1,],y=MatEvol[2,],
         ylim=c(min(MatEvol[-1,]),max(MatEvol[-1,])),
         xlab = info, ylab = "Optimal Criteria",pch=0)
    
    lapply(c(2:nobj),function(i){points(x=MatEvol[1,],y=MatEvol[i+1,],col=i,pch=(i-1))})
    
    xlgd <- popsize +(nrun-popsize)*2/3 ; ylgd <- min(MatEvol[-1,]) + (max(MatEvol[-1,])-min(MatEvol[-1,]))/2
    legend(xlgd,ylgd,legend=objnames,col=1:nobj,pch=0:(nobj-1))
  }
}
