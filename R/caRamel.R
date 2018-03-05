#' caRamel
#'
#' R version of the multi-objective optimizer 'CaRaMEL' originally written for Scilab by Nicolas Le Moine
#'
#' Documentation : "Principe de l'optimiseur CaRaMEL et illustration au travers d'exemples de parametres dans le cadre de la modelisation hydrologique conceptuelle"
#'                 Frederic Hendrickx (EDF) and Nicolas Le Moine (UPMC)
#'                 Report EDF H-P73-2014-09038-FR
#' @param nobj : (integer, length = 1) the number of objectives to optimize (nobj >= 2)
#' @param nvar : (integer, length = 1) the number of variables
#' @param minmax : (logical, length = nobj) the objective is either a minimization (FALSE value) or a maximization (TRUE value)
#' @param bounds : (matrix, nrow = nvar, ncol = 2) lower and upper bounds for the variables
#' @param func : the name of the objective function to optimize. Input argument is the number of parameter set (integer) in the x matrix. The function has to return a vector of at least 'nobj' values (Objectives 1 to nobj are used for optimization, values after nobj are recorded for information.).
#' @param popsize : (integer, length = 1) the population size for the genetic algorithm
#' @param archsize : (integer, length = 1) the size of the Pareto front
#' @param maxrun : (integer, length = 1) the max. number of simulations allowed
#' @param prec : (double, length = nobj) the desired accuracy for the optimization of the objectives
#' @param repart_gene : (integer, length = 4) optional, number of new parameter sets for each rule and per generation
#' @param gpp : (integer, length = 1) optional, calling frequency for the rule "Fireworks"
#' @param blocks (optional): groups for parameters
#' @param pop : (matrix, nrow = nset, ncol = nvar or nvar+nobj ) optional, initial population (used to restart an optimization)
#' @param funcinit (optional): the name of the initialization function applied on each node of cluster when parallel computation. The arguments are cl and numcores.
#' @param noms_obj (optional): the name of the objectives
#' @param listsave (optional): names of the listing files. Default: None (no output). If exists, fields to be defined: "pmt" (file of parameters on the Pareto Front), "obj" (file of corresponding objective values), "evol" (evolution of maximum objectives by generation). Optional field: "totalpop" (total population and corresponding objectives, useful to restart a computation)
#' @param write_gen : (integer, length = 1) optional, if = 1, save files 'pmt' and 'obj' at each generation (= 0 by default)
#' @param carallel : (logical, length = 1) optional, do parallel computations (TRUE by default)
#' @param numcores : (integer, length = 1) optional, the number of cores for the parallel computations (all cores by default).
#
##' @return
##' List of five elements:
##' \describe{
##' \item{success}{return value (logical, length = 1) : TRUE if successfull}
##' \item{parameters}{Pareto front (matrix, nrow = archsize, ncol = nvar)}
##' \item{objectives}{objectives of the Pareto front (matrix, nrow = archsize, ncol = nobj+nadditional)}
##' \item{save_crit}{evolution of the maximum objectives}
##' \item{total_pop}{total population (matrix, nrow = popsize+archsize, ncol = nvar+nobj+nadditional)}
##' }
#' 
#' @examples
#' # Definition of the test function
#' viennet <- function(i) {
#'   val1 <- 0.5*(x[i,1]*x[i,1]+x[i,2]*x[i,2])+sin(x[i,1]*x[i,1]+x[i,2]*x[i,2])
#'   val2 <- 15+(x[i,1]-x[i,2]+1)*(x[i,1]-x[i,2]+1)/27+(3*x[i,1]-2*x[i,2]+4)*(3*x[i,1]-2*x[i,2]+4)/8
#'   val3 <- 1/(x[i,1]*x[i,1]+x[i,2]*x[i,2]+1) -1.1*exp(-(x[i,1]*x[i,1]+x[i,2]*x[i,2]))
#'   return(c(val1,val2,val3))
#' }
#' # Number of objectives
#' nobj <- 3
#' # Number of variables
#' nvar <- 2
#' # All the objectives are to be minimized
#' minmax <- c(FALSE, FALSE, FALSE)
#' # Define the bound constraints
#' bounds <- matrix(data = 1, nrow = nvar, ncol = 2)
#' bounds[, 1] <- -3 * bounds[, 1]
#' bounds[, 2] <- 3 * bounds[, 2]
#' 
#' # Caramel optimization
#' results <-
#'   caRamel(nobj = nobj,
#'           nvar = nvar,
#'           minmax =  minmax,
#'           bounds = bounds,
#'           func = viennet,
#'           popsize = 100,
#'           archsize = 100,
#'           maxrun = 500,
#'           prec = matrix(1.e-3, nrow = 1, ncol = nobj),
#'           carallel = FALSE)
#' 
#' @author Fabrice Zaoui - Celine Monteil

caRamel <-
  function(nobj,
           nvar,
           minmax,
           bounds,
           func,
           popsize,
           archsize,
           maxrun,
           prec,
           repart_gene = c(5, 5, 5, 5),
           gpp = NULL,
           blocks = NULL,
           pop = NULL,
           funcinit = NULL,
           noms_obj = NULL,
           listsave = NULL,
           write_gen = 0,
           carallel = TRUE,
           numcores = NULL) {
    
    # Check the input arguments #
    #############################
    if (nobj <= 1) {
      message("the number of objectives must be greater than one!")
      return(list("success" = FALSE, "message" ="the number of objectives must be greater than one!"))
    }
    if (nvar < 1) {
      message("the number of variables must be greater than zero!")
      return(list("success" = FALSE, "message" ="the number of variables must be greater than zero!"))
    }
    if (length(minmax) != nobj) {
      message("the dimension of 'minmax' is incorrect!")
      return(list("success" = FALSE, "message" ="the dimension of 'minmax' is incorrect!"))
    }
    if ((dim(bounds)[1] != nvar) | (dim(bounds)[2] != 2)) {
      message("the dimension of 'bounds' is incorrect!")
      return(list("success" = FALSE, "message" ="the dimension of 'bounds' is incorrect!"))
    }
    if (class(func) != "function") {
      message("'func' is not an R function!")
      return(list("success" = FALSE, "message" ="'func' is not an R function!"))
    }
    if (popsize < 1) {
      message("'popsize' must be strictly positive!")
      return(list("success" = FALSE, "message" ="'popsize' must be strictly positive!"))
    }
    if (archsize < 1) {
      message("'archsize' must be strictly positive!")
      return(list("success" = FALSE, "message" ="'archsize' must be strictly positive!"))
    }
    if (maxrun < 1) {
      message("'maxrun' must be strictly positive!")
      return(list("success" = FALSE, "message" ="'maxrun' must be strictly positive!"))
    }
    if (length(repart_gene)!=4) {
      message("the dimension of'repart_gene' must be 4!")
      return(list("success" = FALSE, "message" ="the dimension of'repart_gene' must be 4!"))
    }
    if (!is.null(gpp)){
      if (gpp < 1) {
        message("gpp must be greater than zero!")
        return(list("success" = FALSE, "message" ="gpp be greater than zero!"))
      }
    }
    initialise_calc<-0
    if (!is.null(funcinit)){
      if (class(funcinit) != "function") {
        message("'funcinit' is not an R function!")
        return(list("success" = FALSE, "message" ="'funcinit' is not an R function!"))
      }
      initialise_calc<-1
    }
    if (is.null(noms_obj)){noms_obj=paste("Obj",as.character(c(1:nobj)),sep="")}
    writefile<-0
    if (!is.null(listsave)){
      if (class(listsave) != "list") {
        message("'listsave' is not an R list!")
        return(list("success" = FALSE, "message" ="'listsave' is not an R list!"))
      }
      writefile<-1
      if (is.null(listsave$pmt)){
        message(" 'listsave$pmt' must be defined!")
        return(list("success" = FALSE, "message" =" 'listsave$pmt' must be defined!"))
      }
      if (is.null(listsave$obj)){
        message(" 'listsave$obj' must be defined!")
        return(list("success" = FALSE, "message" =" 'listsave$obj' must be defined!"))
      }
      if (is.null(listsave$evol)){
        message(" 'listsave$evol' must be defined!")
        return(list("success" = FALSE, "message" =" 'listsave$evol' must be defined!"))
      }
      ecrit_total_pop = 0 
      if (!is.null(listsave$totalpop)){
        ecrit_total_pop = 1
      }
    }
    if (write_gen==1){
      if (writefile==0){
        message(" 'listsave' must be defined to use write_gen!")
        return(list("success" = FALSE, "message" =" 'listsave' must be defined to use write_gen!"))
      }
      listsave$RadPmt <- gsub(pattern = ".txt",replacement = "",listsave$pmt)
      listsave$RadObj <- gsub(pattern = ".txt",replacement = "",listsave$obj)
      if (ecrit_total_pop==1){listsave$RadPop <- gsub(pattern = ".txt",replacement = "",listsave$totalpop)}
    }
    if (typeof(carallel)!="logical") {
      message("'carallel' must be a logical!")
      return(list("success" = FALSE, "message" ="'carallel' must be a logical!"))
    }
    
    # Initializations
    save_crit<<-c()
    sp <- (bounds[, 2] - bounds[, 1]) / (2 * sqrt(3)) # standard deviation
    gsearch <-
      ceiling((nvar / 10) * nobj / log(nobj + 1)) # independant search every 'gsearch' iteration
    precis <-
      matrix(data = prec, nrow = nobj) # precision tolerance for each objective
    nrun <- 0 # main loop index
    ngen <- 0 # generation number
    if (is.null(gpp)){
      gpp <- ceiling( nvar * (nobj+1) * 4 / sum(repart_gene) )
    }
    
    # Init the parallel computation
    if (carallel==TRUE){
      if (is.null(numcores)){ numcores <- detectCores()}
      cl <- makeCluster(numcores)
      if (initialise_calc==1){
        funcinit(cl,numcores) # Init for "func"
      }
    }
    
    # Check the type of the initial population
    if (!is.null(pop)){
      pop <- as.matrix(pop)
      if(length(pop[1,])<(nvar+nobj)){  # If an initial population exists with no corresponding values for the objectives
        
        # Evaluation of the objectives
        x<<-pop[,1:nvar]
        if (carallel==TRUE){
          newfeval <- NULL
          clusterExport(cl=cl, varlist=c("x"), envir = environment())
          res = parLapply(cl, 1:dim(x)[1], func)
          for (j in 1:dim(x)[1]) {
            newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
          }
        } else { #sequential calls
          newfeval <- matrix(data = 0.,
                             nrow = dim(x)[1],
                             ncol = nobj)
          for (i in 1:dim(x)[1]) {
            res <- func(i)
            newfeval[i, ] <- res[1:nobj]
          }
        }
        pop <- cbind(pop[,1:nvar],newfeval)
        nrun <- nrun + length(pop[,1])
      }
    }
    
    # Optimization
    message("Beginning of optimization process")
    pb <- txtProgressBar(min=0, max=1,initial=0,title="caRamel progress :",label="caRamel progress :" , style=3)
    while (nrun < maxrun) {
      ngen <- ngen + 1
      
      # new population
      if (is.null(pop)) {
        A <-
          t(matrix(data = bounds[, 2] - bounds[, 1],
                   ncol = popsize,
                   nrow = nvar))
        B <- t(matrix(data = bounds[, 1],
                      ncol = popsize,
                      nrow = nvar))
        x <- A * matrix(runif(popsize * nvar), ncol = nvar) + B
        probj <- matrix(data = NaN,
                        nrow = popsize,
                        ncol = nobj)
      } else {
        vamax <- (ngen %% gpp) == 0
        Xp <-
          newXval(as.matrix(pop[, 1:nvar]),
                  pop[, (nvar + 1):(nvar + nobj)],
                  minmax,
                  sp,
                  bounds,
                  repart_gene,
                  blocks,
                  vamax)
        x <- Xp$x
        probj <- Xp$pcrit
      }
      
      # simulations
      additional_eval <- NULL
      # parallel calls
      if (carallel==TRUE){
        newfeval <- NULL
        clusterExport(cl=cl, varlist=c("x"), envir = environment())
        res <- parLapply(cl, 1:dim(x)[1], func)
        for (j in 1:dim(x)[1]) {
          newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
          additional_eval <- rbind(additional_eval, res[[j]][(nobj+1):length(res[[1]])])
        }
        nadditional <- ncol(additional_eval)
        if (length(res[[1]])<(nobj+1)){
          additional_eval <- NULL
          nadditional <- 0
        }
      } else { #sequential calls
        newfeval <- matrix(data = 0.,
                           nrow = dim(x)[1],
                           ncol = nobj)
        x<<-x
        for (i in 1:dim(x)[1]) {
          res <- func(i)
          newfeval[i, ] <- res[1:nobj]
          additional_eval <- rbind(additional_eval, res[(nobj+1):length(res)])
        }
        nadditional <- ncol(additional_eval)
        if (length(res)<(nobj+1)){
          additional_eval <- NULL
          nadditional <- 0
        }
      }
      
      
      # update the number of calls
      nrun <- nrun + dim(x)[1]
      
      # show progress
      progress <- round(min(nrun,maxrun) * 100 / maxrun)
      setTxtProgressBar(pb, min(nrun,maxrun)/maxrun)
      
      # deal with NaN value
      detect_nan <- is.na(newfeval)
      set_ok <- !rowSums(detect_nan)
      newfeval <- newfeval[set_ok, ]
      x <- x[set_ok, ]
      probj <- probj[set_ok, ]
      pop1 <- rbind(pop, cbind(x, newfeval, additional_eval))
      
      # decrease population
      #message("decrease pop")
      ind <-
        decrease_pop(pop1[, (nvar + 1):(nvar + nobj)], minmax, prec, archsize, popsize)
      
      # archive
      # message("archive")
      arch <- matrix(pop1[ind$arch, ], nrow=length(ind$arch), ncol=nobj+nvar+nadditional)
      
      # population update
      #message("pop update")
      pop <- pop1[c(ind$arch, ind$pop), ]
      param_arch <- arch[, 1:nvar]
      crit_arch <- matrix(arch[, (nvar + 1):(nvar + nobj)], nrow=length(ind$arch), ncol=nobj)
      if (nadditional>0){
        additional_eval[,1:ncol(additional_eval)] <- additional_eval[set_ok, 1:nadditional]
        additional_eval <- matrix(arch[, (nvar + nobj+1):(nvar + nobj+nadditional)], nrow=length(ind$arch), ncol=nadditional)
      }
      
      # Records on criteria
      a=c(lapply(c(1:nobj),function(i){max(crit_arch[,i])}))
      maxcrit=as.data.frame(a,col.names = noms_obj)
      a=c(lapply(c(1:nobj),function(i){min(crit_arch[,i])}))
      mincrit=as.data.frame(a,col.names = noms_obj)
      crit=mincrit; crit[minmax]<-maxcrit[minmax]
      save_crit<-cbind(save_crit,c(nrun,t(crit)))
      
      # Graphs
      info<-paste("ngen=",ngen,", nrun=",nrun,", gpp=",gpp,sep="")
      nbre_fen = choose(n = nobj, k=2)+1
      if (nbre_fen <= 4){
        par(mfrow=c(2,2))
      } else {
        par(mfrow = c(3,floor(nbre_fen/3)+1))
      }
      l = seq(1,nobj)
      for (i_fig in 1:(nobj-1)){
        l_tmp = l[-i_fig]; l_tmp = l_tmp[l_tmp>i_fig]
        for (i_fig2 in l_tmp){
          plot(crit_arch[, i_fig], crit_arch[, i_fig2],xlab = noms_obj[i_fig],ylab = noms_obj[i_fig2])
        }
      }
      plot(x=save_crit[1,],y=save_crit[2,], ylim=c(min(save_crit[-1,]),max(save_crit[-1,])),xlab = info, ylab = "Optimal Criteria")
      lapply(c(2:nobj),function(i){points(x=save_crit[1,],y=save_crit[i+1,],col=i)})
      xlgd <- popsize +(nrun-popsize)*2/3 ; ylgd <- min(save_crit[-1,]) + (max(save_crit[-1,])-min(save_crit[-1,]))/2
      legend(xlgd,ylgd,legend=noms_obj,col=1:nobj,fill=1:nobj)
      
      # Online saves
      if (writefile == 1){
        if (write_gen == 1){
          listsave$pmt <- paste(listsave$RadPmt,"_gen",ngen,".txt",sep="")
          listsave$obj <- paste(listsave$RadObj,"_gen",ngen,".txt",sep="")
          if (ecrit_total_pop==1){
            listsave$totalpop <- paste(listsave$RadPop,"_gen",ngen,".txt",sep="")
          }
        }
        
        write.table(param_arch,listsave$pmt,row.names = FALSE,col.names = FALSE)
        write.table(cbind(crit_arch,additional_eval),listsave$obj,row.names = FALSE,col.names = FALSE)
        write.table(t(save_crit),listsave$evol,row.names = FALSE,col.names = FALSE)
        if (ecrit_total_pop==1){
          write.table(pop,listsave$totalpop,row.names = FALSE,col.names = FALSE)
        }
      }
    }
    
    if (carallel==TRUE){stopCluster(cl)}
    close(pb)
    
    return(list(
      "success" = TRUE,
      "parameters" = param_arch,
      "objectives" = cbind(crit_arch,additional_eval),
      "save_crit" = t(save_crit),
      "total_pop"= pop
    ))
  }
