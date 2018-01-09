#' caRamel
#'
#' R version of the multi-objective optimizer 'CaRaMEL' originally written for Scilab by Nicolas Le Moine
#'
#' Documentation : "Principe de l'optimiseur CaRaMEL et illustration au travers d'exemples de parametres dans le cadre de la modelisation hydrologique conceptuelle"
#'                 Frederic Hendrickx (EDF) and Nicolas Le Moine (UPMC)
#'                 Report EDF H-P73-2014-09038-FR
#' @param nobj (integer, length = 1) the number of objectives to optimize (nobj >= 2)
#' @param nvar (integer, length = 1) the number of variables
#' @param minmax (logical, length = nobj) the objective is either a minimization (FALSE value) or a maximization (TRUE value)
#' @param bounds (matrix, nrow = nvar, ncol = 2) lower and upper bounds for the variables
#' @param func the name of the objective function to optimize. The function has to return 'nobj' values with : nobj >= 2
#' @param popsize (integer, length = 1) the population size for the genetic algorithm
#' @param archsize (integer, length = 1) the size of the Pareto front
#' @param maxrun (integer, length = 1) the max. number of simulations allowed
#' @param prec (double, length = nobj) the desired precision for the optimization of the objectives
#' @param repart_gene (integer, length = 4) nombre de jeux generes pour chaque regle par generation,
#' @param gpp (integer, length = 1) frequence d'appel de la regle "Fireworks"
#' @param blocs blocs pour les parametres
#' @param pop (matrix, nrow = njeux, ncol = nvar ou nvar+nobj ) Population initiale
#' @param funcinit the name of the initialization function. The arguments are is cl and numcores.
#' @param noms_obj the name of the objectives
#' @param listsave liste du noms des fichiers de suivi ecrits en cours d'optim, pas d'ecriture par defaut. Champs obligatoires : "pmt" (fichier des jeux de parametres sur le front), "obj" (Objectifs associes), "evol" (evolution des objectifs max par generation). Champ optionnel : "totalpop" (population totale : parametres et objectifs associes)
#' @param ecrit_gen (integer, length = 1) si =1, conservation des fichiers pmt et obj a chaque generation (0 par defaut)
#' @param parallele (logical, length = 1) activation du calcul parallele (TRUE par defaut)
#' @param numcores (integer, length = 1) the number of cores for parallel computation
#
##' @return
##' List of five elements:
##' \describe{
##' \item{success}{return value (logical, length = 1) : TRUE if successfull}
##' \item{parameters}{Pareto front (matrix, nrow = archsize, ncol = nvar)}
##' \item{objectives}{objectives of the Pareto front (matrix, nrow = archsize, ncol = nobj)}
##' \item{suivi_crit}{evolution of the maximum objectives}
##' \item{pop_totale}{total population (matrix, nrow = popsize+archsize, ncol = nvar+nobj+ncomplement)}
##' }
#' @author Fabrice Zaoui - Celine Monteil
#' @export

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
           blocs = NULL,
           pop = NULL,
           funcinit = NULL,
           noms_obj = NULL,
           listsave = NULL,
           ecrit_gen = 0,
           parallele = TRUE,
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
    ecritfic<-0
    if (!is.null(listsave)){
      if (class(listsave) != "list") {
        message("'listsave' is not an R list!")
        return(list(-1, "'listsave' is not an R list!", NA, NA))
      }
      ecritfic<-1
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
      ecrit_pop_totale = 0 
      if (!is.null(listsave$totalpop)){
        ecrit_pop_totale = 1
      }
    }
    if (ecrit_gen==1){
      if (ecritfic==0){
        message(" 'listsave' must be defined to use ecrit_gen!")
        return(list(-1, " 'listsave' must be defined to use ecrit_gen!", NA, NA))
      }
      listsave$RadPmt <- gsub(pattern = ".txt",replacement = "",listsave$pmt)
      listsave$RadObj <- gsub(pattern = ".txt",replacement = "",listsave$obj)
      if (ecrit_pop_totale==1){listsave$RadPop <- gsub(pattern = ".txt",replacement = "",listsave$totalpop)}
    }
    if (typeof(parallele)!="logical") {
      message("'parallele' must be a logical!")
      return(list(-1, "'parallele' must be a logical!", NA, NA))
    }
    
    # Initializations
    suivi_crit<<-c()
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
    
    # Initialisation du calcul parallele
    if (parallele==TRUE){
      if (is.null(numcores)){ numcores <- detectCores()}
      cl <- makeCluster(numcores)
      if (initialise_calc==1){
        funcinit(cl,numcores) # Initialisation pour le simulateur "func"
      }
    }
    
    # Verification du type de la population initiale
    if (!is.null(pop)){
      pop <- as.matrix(pop)
      if(length(pop[1,])<(nvar+nobj)){  # S'il y a une pop initiale mais que les objectifs ne sont pas evalues
        
        # Evaluation des objectifs
        x<<-pop[,1:nvar]
        eval_complementaire <- NULL
        if (parallele==TRUE){
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
                  blocs,
                  vamax)
        x <- Xp$x
        probj <- Xp$pcrit
      }
      
      # simulations
      eval_complementaire <- NULL
      # parallel calls
      if (parallele==TRUE){
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
      
      # Suivi criteres
      a=c(lapply(c(1:nobj),function(i){max(crit_arch[,i])}))
      maxcrit=as.data.frame(a,col.names = noms_obj)
      a=c(lapply(c(1:nobj),function(i){min(crit_arch[,i])}))
      mincrit=as.data.frame(a,col.names = noms_obj)
      crit=mincrit; crit[minmax]<-maxcrit[minmax]
      suivi_crit<-cbind(suivi_crit,c(nrun,t(crit)))
      
      # Graphes
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
      plot(x=suivi_crit[1,],y=suivi_crit[2,], ylim=c(min(suivi_crit[-1,]),max(suivi_crit[-1,])),xlab = info, ylab = "Criteres Optimaux")
      lapply(c(2:nobj),function(i){points(x=suivi_crit[1,],y=suivi_crit[i+1,],col=i)})
      xlgd <- popsize +(nrun-popsize)*2/3 ; ylgd <- min(suivi_crit[-1,]) + (max(suivi_crit[-1,])-min(suivi_crit[-1,]))/2
      legend(xlgd,ylgd,legend=noms_obj,col=1:nobj,fill=1:nobj)
      
      # Pour sauvegarde en cours d'optimisation
      if (ecritfic == 1){
        if (ecrit_gen == 1){
          listsave$pmt <- paste(listsave$RadPmt,"_gen",ngen,".txt",sep="")
          listsave$obj <- paste(listsave$RadObj,"_gen",ngen,".txt",sep="")
          if (ecrit_pop_totale==1){
            listsave$totalpop <- paste(listsave$RadPop,"_gen",ngen,".txt",sep="")
          }
        }
        
        write.table(param_arch,listsave$pmt,row.names = FALSE,col.names = FALSE)
        write.table(cbind(crit_arch,eval_complementaire),listsave$obj,row.names = FALSE,col.names = FALSE)
        write.table(t(suivi_crit),listsave$evol,row.names = FALSE,col.names = FALSE)
        if (ecrit_pop_totale==1){
          write.table(pop,listsave$totalpop,row.names = FALSE,col.names = FALSE)
        }
      }
    }
    
    if (parallele==TRUE){stopCluster(cl)}
    close(pb)
    
    return(list(
      "success" = TRUE,
      "parameters" = param_arch,
      "objectives" = cbind(crit_arch,eval_complementaire),
      "suivi_crit" = t(suivi_crit),
      "pop_totale"= pop
    ))
  }
