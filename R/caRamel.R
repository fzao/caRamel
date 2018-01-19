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
#' @param func : the name of the objective function to optimize. The function has to return 'nobj' values with : nobj >= 2
#' @param popsize : (integer, length = 1) the population size for the genetic algorithm
#' @param archsize : (integer, length = 1) the size of the Pareto front
#' @param maxrun : (integer, length = 1) the max. number of simulations allowed
#' @param prec : (double, length = nobj) the desired accuracy for the optimization of the objectives
#' @param repart_gene : (integer, length = 4) number of new solutions for each rule and per generation
#' @param gpp : (integer, length = 1) calling frequency for the rule "Fireworks"
#' @param blocks : groups for parameters
#' @param pop : (matrix, nrow = nset, ncol = nvar or nvar+nobj ) Initial population
#' @param funcinit : the name of the initialization function. The arguments are cl and numcores
#' @param noms_obj : the name of the objectives
#' @param listsave : names of the listing files. Default: None (no output)
#' @param write_gen : (integer, length = 1) if = 1, save files 'pmt' and 'obj' at each generation (= 0 by default)
#' @param carallel : (logical, length = 1) do parallel computations (TRUE by default)
#' @param numcores : (integer, length = 1) the number of cores for the parallel computations
#
##' @return
##' List of five elements:
##' \describe{
##' \item{success}{return value (logical, length = 1) : TRUE if successfull}
##' \item{parameters}{Pareto front (matrix, nrow = archsize, ncol = nvar)}
##' \item{objectives}{objectives of the Pareto front (matrix, nrow = archsize, ncol = nobj)}
##' \item{save_crit}{evolution of the maximum objectives}
##' \item{total_pop}{total population (matrix, nrow = popsize+archsize, ncol = nvar+nobj+ncomplement)}
##' }
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
      return(list(-1, "the number of objectives must be greater than one!", NA, NA))
    }
    if (nvar < 1) {
      message("the number of variables must be greater than zero!")
      return(list(-1, "the number of variables must be greater than zero!", NA, NA))
    }
    if (length(minmax) != nobj) {
      message("the dimension of 'minmax' is incorrect!")
      return(list(-1, "the dimension of 'minmax' is incorrect!", NA, NA))
    }
    if ((dim(bounds)[1] != nvar) | (dim(bounds)[2] != 2)) {
      message("the dimension of 'bounds' is incorrect!")
      return(list(-1, "the dimension of 'bounds' is incorrect!", NA, NA))
    }
    if (class(func) != "function") {
      message("'func' is not an R function!")
      return(list(-1, "'func' is not an R function!", NA, NA))
    }
    if (popsize < 1) {
      message("'popsize' must be strictly positive!")
      return(list(-1, "'popsize' must be strictly positive!", NA, NA))
    }
    if (archsize < 1) {
      message("'archsize' must be strictly positive!")
      return(list(-1, "'archsize' must be strictly positive!", NA, NA))
    }
    if (maxrun < 1) {
      message("'maxrun' must be strictly positive!")
      return(list(-1, "'maxrun' must be strictly positive!", NA, NA))
    }
    if (length(repart_gene)!=4) {
      message("the dimension of'repart_gene' must be 4!")
      return(list(-1, "the dimension of'repart_gene' must be 4!", NA, NA))
    }
    if (!is.null(gpp)){
      if (gpp < 1) {
        message("gpp must be greater than zero!")
        return(list(-1, "gpp be greater than zero!", NA, NA))
      }
    }
    initialise_calc<-0
    if (!is.null(funcinit)){
      if (class(funcinit) != "function") {
        message("'funcinit' is not an R function!")
        return(list(-1, "'funcinit' is not an R function!", NA, NA))
      }
      initialise_calc<-1
    }
    if (is.null(noms_obj)){noms_obj=paste("Obj",as.character(c(1:nobj)),sep="")}
    writefile<-0
    if (!is.null(listsave)){
      if (class(listsave) != "list") {
        message("'listsave' is not an R list!")
        return(list(-1, "'listsave' is not an R list!", NA, NA))
      }
      writefile<-1
      if (is.null(listsave$pmt)){
        message(" 'listsave$pmt' must be defined!")
        return(list(-1, " 'listsave$pmt' must be defined!", NA, NA))
      }
      if (is.null(listsave$obj)){
        message(" 'listsave$obj' must be defined!")
        return(list(-1, " 'listsave$obj' must be defined!", NA, NA))
      }
      if (is.null(listsave$evol)){
        message(" 'listsave$evol' must be defined!")
        return(list(-1, " 'listsave$evol' must be defined!", NA, NA))
      }
      ecrit_total_pop = 0 
      if (!is.null(listsave$totalpop)){
        ecrit_total_pop = 1
      }
    }
    if (write_gen==1){
      if (writefile==0){
        message(" 'listsave' must be defined to use write_gen!")
        return(list(-1, " 'listsave' must be defined to use write_gen!", NA, NA))
      }
      listsave$RadPmt <- gsub(pattern = ".txt",replacement = "",listsave$pmt)
      listsave$RadObj <- gsub(pattern = ".txt",replacement = "",listsave$obj)
      if (ecrit_total_pop==1){listsave$RadPop <- gsub(pattern = ".txt",replacement = "",listsave$totalpop)}
    }
    if (typeof(carallel)!="logical") {
      message("'carallel' must be a logical!")
      return(list(-1, "'carallel' must be a logical!", NA, NA))
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
        eval_complementaire <- NULL
        if (carallel==TRUE){
          newfeval <- NULL
          clusterExport(cl=cl, varlist=c("x"), envir = environment())
          res = parLapply(cl, 1:dim(x)[1], func)
          for (j in 1:dim(x)[1]) {
            newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
            eval_complementaire <- rbind(eval_complementaire, res[[j]][(nobj+1):length(res[[1]])])
          }
        } else { #sequential calls
          newfeval <- matrix(data = 0.,
                             nrow = dim(x)[1],
                             ncol = nobj)
          for (i in 1:dim(x)[1]) {
            res <- func(i)
            newfeval[i, ] <- res[1:nobj]
            eval_complementaire <- rbind(eval_complementaire, res[(nobj+1):length(res)])
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
      eval_complementaire <- NULL
      # parallel calls
      if (carallel==TRUE){
        newfeval <- NULL
        clusterExport(cl=cl, varlist=c("x"), envir = environment())
        res <- parLapply(cl, 1:dim(x)[1], func)
        for (j in 1:dim(x)[1]) {
          newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
          eval_complementaire <- rbind(eval_complementaire, res[[j]][(nobj+1):length(res[[1]])])
        }
        ncomplement <- ncol(eval_complementaire)
        if (length(res[[1]])<(nobj+1)){
          eval_complementaire <- NULL
          ncomplement <- 0
        }
      } else { #sequential calls
        newfeval <- matrix(data = 0.,
                           nrow = dim(x)[1],
                           ncol = nobj)
        x<<-x
        for (i in 1:dim(x)[1]) {
          res <- func(i)
          newfeval[i, ] <- res[1:nobj]
          eval_complementaire <- rbind(eval_complementaire, res[(nobj+1):length(res)])
        }
        ncomplement <- ncol(eval_complementaire)
        if (length(res)<(nobj+1)){
          eval_complementaire <- NULL
          ncomplement <- 0
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
      pop1 <- rbind(pop, cbind(x, newfeval, eval_complementaire))
      
      # decrease population
      #message("decrease pop")
      ind <-
        decrease_pop(pop1[, (nvar + 1):(nvar + nobj)], minmax, prec, archsize, popsize)
      
      # archive
      # message("archive")
      arch <- matrix(pop1[ind$arch, ], nrow=length(ind$arch), ncol=nobj+nvar+ncomplement)
      
      # population update
      #message("pop update")
      pop <- pop1[c(ind$arch, ind$pop), ]
      param_arch <- arch[, 1:nvar]
      crit_arch <- matrix(arch[, (nvar + 1):(nvar + nobj)], nrow=length(ind$arch), ncol=nobj)
      if (ncomplement>0){
        eval_complementaire[,1:ncol(eval_complementaire)] <- eval_complementaire[set_ok, 1:ncomplement]
        eval_complementaire <- matrix(arch[, (nvar + nobj+1):(nvar + nobj+ncomplement)], nrow=length(ind$arch), ncol=ncomplement)
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
        write.table(cbind(crit_arch,eval_complementaire),listsave$obj,row.names = FALSE,col.names = FALSE)
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
      "objectives" = cbind(crit_arch,eval_complementaire),
      "save_crit" = t(save_crit),
      "total_pop"= pop
    ))
  }
