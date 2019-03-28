#' Plotting of caRamel results
#' 
#' Plot graphs of the Pareto front and a graph of optimization evolution
#' 
#' @param caramel_results : list resulting from the caRamel() function, with fields $objectives and $save_crit
#' @param nobj : number of objectives (optional)
#' @param objnames : vector of objectives names (optional)
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
#' plot_caramel(results)

plot_caramel <- function(caramel_results, nobj=NULL, objnames=NULL){
  
  ngen <- length(caramel_results$save_crit[,1])
  nrun <- caramel_results$save_crit[ngen,1]
  
  if (is.null(nobj)){nobj <- length(caramel_results$objectives[1,])}
  if (is.null(objnames)){objnames <- paste("Obj",as.character(c(1:nobj)),sep="")}else{nobj<-length(objnames)}
  
  MatObj <- caramel_results$objectives
  MatEvol <- t(caramel_results$save_crit)
  
  plot_population(MatObj, nobj, ngen, nrun, objnames, MatEvol)
  
}
