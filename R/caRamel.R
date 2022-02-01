#' MAIN FUNCTION: multi-objective optimizer
#'
#' Multi-objective optimizer. It requires to define a multi-objective function (func) to calibrate the model and bounds on the parameters to optimize.
#'
#' The optimizer was originally written for Scilab by Nicolas Le Moine.
#' The algorithm is a hybrid of the MEAS algorithm (Efstratiadis and Koutsoyiannis (2005) <doi:10.13140/RG.2.2.32963.81446>) by using the directional search method based on the simplexes of the objective space
#'     and the epsilon-NGSA-II algorithm with the method of classification of the parameter vectors archiving management by epsilon-dominance (Reed and Devireddy <doi:10.1142/9789812567796_0004>).
#' Reference : "Multi-objective calibration by combination of stochastic and gradient-like parameter generation rules – the caRamel algorithm"
#'             Celine Monteil (EDF), Fabrice Zaoui (EDF), Nicolas Le Moine (UPMC) and Frederic Hendrickx (EDF)
#'             June 2020 Hydrology and Earth System Sciences 24(6):3189-3209
#'             DOI: 10.5194/hess-24-3189-2020
#' Documentation : "Principe de l'optimiseur CaRaMEL et illustration au travers d'exemples de parametres dans le cadre de la modelisation hydrologique conceptuelle"
#'                 Frederic Hendrickx (EDF) and Nicolas Le Moine (UPMC)
#'                 Report EDF H-P73-2014-09038-FR
#'
#' @param nobj : (integer, length = 1) the number of objectives to optimize (nobj >= 2)
#' @param nvar : (integer, length = 1) the number of variables
#' @param minmax : (logical, length = nobj) the objective is either a minimization (FALSE value) or a maximization (TRUE value)
#' @param bounds : (matrix, nrow = nvar, ncol = 2) lower and upper bounds for the variables
#' @param func : (function) the objective function to optimize. Input argument is the number of parameter set (integer) in the x matrix. The function has to return a vector of at least 'nobj' values (Objectives 1 to nobj are used for optimization, values after nobj are recorded for information.).
#' @param popsize : (integer, length = 1) the population size for the genetic algorithm
#' @param archsize : (integer, length = 1) the size of the Pareto front
#' @param maxrun : (integer, length = 1) the max. number of simulations allowed
#' @param prec : (double, length = nobj) the desired accuracy for the optimization of the objectives
#' @param repart_gene : (integer, length = 4) optional, number of new parameter sets for each rule and per generation
#' @param gpp : (integer, length = 1) optional, calling frequency for the rule "Fireworks"
#' @param blocks (optional): groups for parameters
#' @param pop : (matrix, nrow = nset, ncol = nvar or nvar+nobj ) optional, initial population (used to restart an optimization)
#' @param funcinit (function, optional): the initialization function applied on each node of cluster when parallel computation. The arguments are cl and numcores
#' @param objnames (optional): names of the objectives
#' @param listsave (optional): names of the listing files. Default: None (no output). If exists, fields to be defined: "pmt" (file of parameters on the Pareto Front), "obj" (file of corresponding objective values), "evol" (evolution of maximum objectives by generation). Optional field: "totalpop" (total population and corresponding objectives, useful to restart a computation)
#' @param write_gen : (logical, length = 1) optional, if TRUE, save files 'pmt' and 'obj' at each generation (FALSE by default)
#' @param carallel : (integer, length = 1) optional, do parallel computations? (0: sequential, 1:parallel (default) , 2:user-defined choice)
#' @param numcores : (integer, length = 1) optional, the number of cores for the parallel computations (all cores by default)
#' @param graph : (logical, length = 1) optional, plot graphical output at each generation (TRUE by default)
#' @param sensitivity : (logical, length = 1) optional, compute the first order derivatives of the pareto front (FALSE by default)
#' @param verbose : (logical, length = 1) optional, verbosity mode (TRUE by default)
#' @param worklist : optional values to be transmitted to the user's function (not used)
#
##' @return
##' List of seven elements:
##' \describe{
##' \item{success}{return value (logical, length = 1) : TRUE if successfull}
##' \item{parameters}{Pareto front (matrix, nrow = archsize, ncol = nvar)}
##' \item{objectives}{objectives of the Pareto front (matrix, nrow = archsize, ncol = nobj+nadditional)}
##' \item{derivatives}{list of the Jacobian matrices of the Pareto front if the sensitivity parameter is TRUE or NA otherwise}
##' \item{save_crit}{evolution of the optimal objectives}
##' \item{total_pop}{total population (matrix, nrow = popsize+archsize, ncol = nvar+nobj+nadditional)}
##' \item{gpp}{the calling period for the third generation rule (independent sampling with a priori parameters variance)}
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
#'           carallel = 0)
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
           objnames = NULL,
           listsave = NULL,
           write_gen = FALSE,
           carallel = 1,
           numcores = NULL,
           graph = TRUE,
           sensitivity = FALSE,
           verbose = TRUE,
           worklist = NULL) {

    start_time <- Sys.time()

    # Check the input arguments #
    #############################
    if (nobj <= 1) {
      if (verbose) message("the number of objectives must be greater than one!")
      return(list("success" = FALSE, "message" = "the number of objectives must be greater than one!"))
    }
    if (nvar < 1) {
      if (verbose) message("the number of variables must be greater than zero!")
      return(list("success" = FALSE, "message" = "the number of variables must be greater than zero!"))
    }
    if (length(minmax) != nobj) {
      if (verbose) message("the dimension of 'minmax' is incorrect!")
      return(list("success" = FALSE, "message" = "the dimension of 'minmax' is incorrect!"))
    }
    if ((dim(bounds)[1] != nvar) | (dim(bounds)[2] != 2)) {
      if (verbose) message("the dimension of 'bounds' is incorrect!")
      return(list("success" = FALSE, "message" = "the dimension of 'bounds' is incorrect!"))
    }
    if (class(func) != "function") {
      if (verbose) message("'func' is not an R function!")
      return(list("success" = FALSE, "message" = "'func' is not an R function!"))
    }
    if (popsize < 1) {
      if (verbose) message("'popsize' must be strictly positive!")
      return(list("success" = FALSE, "message" = "'popsize' must be strictly positive!"))
    }
    if (archsize < 1) {
      if (verbose) message("'archsize' must be strictly positive!")
      return(list("success" = FALSE, "message" = "'archsize' must be strictly positive!"))
    }
    if (maxrun < 1) {
      if (verbose) message("'maxrun' must be strictly positive!")
      return(list("success" = FALSE, "message" = "'maxrun' must be strictly positive!"))
    }
    if (length(repart_gene) != 4) {
      if (verbose) message("the dimension of'repart_gene' must be 4!")
      return(list("success" = FALSE, "message" = "the dimension of 'repart_gene' must be 4!"))
    }
    if (sum(repart_gene <= 0) > 0) {
      if (verbose) message("parameter values for each rule of 'repart_gene' must be strictly positive!")
      return(list("success" = FALSE, "message" = "parameter values for each rule of 'repart_gene' must be strictly positive!"))
    }
    if (!is.null(gpp)) {
      if (gpp < 1) {
        if (verbose) message("gpp must be greater than zero!")
        return(list("success" = FALSE, "message" = "gpp be greater than zero!"))
      }
    }
    initialise_calc <- 0
    if (!is.null(funcinit)) {
      if (class(funcinit) != "function") {
        if (verbose) message("'funcinit' is not an R function!")
        return(list("success" = FALSE, "message" = "'funcinit' is not an R function!"))
      }
      initialise_calc <- 1
    }
    if (is.null(objnames)) objnames <- paste("Obj", as.character(c(1:nobj)), sep = "")
    writefile <- FALSE
    if (!is.null(listsave)) {
      if (class(listsave) != "list") {
        if (verbose) message("'listsave' is not an R list!")
        return(list("success" = FALSE, "message" = "'listsave' is not an R list!"))
      }
      writefile <- TRUE
      if (is.null(listsave$pmt)) {
        if (verbose) message(" 'listsave$pmt' must be defined!")
        return(list("success" = FALSE, "message" = " 'listsave$pmt' must be defined!"))
      }
      if (is.null(listsave$obj)) {
        if (verbose) message(" 'listsave$obj' must be defined!")
        return(list("success" = FALSE, "message" = " 'listsave$obj' must be defined!"))
      }
      if (is.null(listsave$evol)) {
        if (verbose) message(" 'listsave$evol' must be defined!")
        return(list("success" = FALSE, "message" = " 'listsave$evol' must be defined!"))
      }
      ecrit_total_pop <- 0
      if (!is.null(listsave$totalpop)) {
        ecrit_total_pop <- 1
      }
    }
    if (write_gen == TRUE) {
      if (writefile == FALSE) {
        if (verbose) message(" 'listsave' must be defined to use write_gen!")
        return(list("success" = FALSE, "message" = " 'listsave' must be defined to use write_gen!"))
      }
      listsave$RadPmt <- gsub(pattern = ".txt", replacement = "", listsave$pmt)
      listsave$RadObj <- gsub(pattern = ".txt", replacement = "", listsave$obj)
      if (ecrit_total_pop == 1) listsave$RadPop <- gsub(pattern = ".txt", replacement = "", listsave$totalpop)
    }
    carallel <- as.integer(carallel)
    if (!(carallel %in% c(0, 1, 2))) {
      if (verbose) message("'carallel' value is not 0, 1 or 2")
      return(list("success" = FALSE, "message" = "wrong value for paramemter 'carallel'"))
    }

    # Initializations
    save_crit <<- c()
    sp <- (bounds[, 2] - bounds[, 1]) / (2 * sqrt(3)) # standard deviation
    gsearch <- ceiling((nvar / 10) * nobj / log(nobj + 1)) # independant search every 'gsearch' iteration
    precis <- matrix(data = prec, nrow = nobj) # precision tolerance for each objective
    nrun <- 0 # main loop index
    ngen <- 0 # generation number
    if (is.null(gpp)) gpp <- ceiling(nvar * (nobj + 1) * 4 / sum(repart_gene))
    jacobian <- NA

    # Init the parallel computation
    if (carallel == 1) {
      if (is.null(numcores)) numcores <- detectCores()
      cl <- makeCluster(numcores)
      if (initialise_calc == 1) funcinit(cl, numcores) # Init for "func"
    }

    # Check the type of the initial population
    if (!is.null(pop)) {
      pop <- as.matrix(pop)
      if(length(pop[1, ]) < (nvar + nobj)) {  # If an initial population exists with no corresponding values for the objectives

        # Evaluation of the objectives
        x <<- pop[, 1:nvar]
        if (carallel == 1) {
          newfeval <- NULL
          clusterExport(cl = cl, varlist = c("x", "worklist"), envir = environment())
          res <- parLapply(cl, 1:dim(x)[1], func)
          for (j in 1:dim(x)[1]) {
            newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
          }
        } else if (carallel == 0) { #sequential calls
          newfeval <- matrix(data = 0.,
                             nrow = dim(x)[1],
                             ncol = nobj)
          for (i in 1:dim(x)[1]) {
            res <- func(i)
            newfeval[i, ] <- res[1:nobj]
          }
        } else {
          newfeval <- func(pop[, 1:nvar])
        }
        pop <- cbind(pop[, 1:nvar], newfeval)
        nrun <- nrun + length(pop[, 1])
      }
    }

    # Optimization
    if(verbose) {
      message(paste("Beginning of caRamel optimization <--", date()))
      message(paste("Number of variables :", as.character(nvar)))
      message(paste("Number of functions :", as.character(nobj)))
      pb <- txtProgressBar(min = 0, max = 1, initial = 0, title = "caRamel progress :", label = "caRamel progress :" , style = 3)
    }
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
        param <- as.matrix(pop[, 1:nvar])
        dim(param) <- c(dim(pop)[1], nvar)
        if (dim(param)[1] < 4) {
          if (carallel == 1) stopCluster(cl)
          if (verbose) {
            close(pb)
            message("Optimization failed")
          }
          return(list("success" = FALSE, "message" = "The number of feasible points is not sufficient! Try to increase the size of the population..."))
        }
        crit <- pop[, (nvar + 1):(nvar + nobj)]
        dim(crit) <- c(dim(pop)[1], nobj)
        Xp <-
          newXval(param,
                  crit,
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
      if (carallel == 1) {
        newfeval <- NULL
        clusterExport(cl = cl, varlist = c("x", "worklist"), envir = environment())
        res <- parLapply(cl, 1:dim(x)[1], func)
        for (j in 1:dim(x)[1]) {
          newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
          additional_eval <- rbind(additional_eval, res[[j]][(nobj+1):length(res[[1]])])
        }
        nadditional <- ncol(additional_eval)
        if (length(res[[1]]) < (nobj + 1)) {
          additional_eval <- NULL
          nadditional <- 0
        }
      } else if (carallel == 0) { #sequential calls
        newfeval <- matrix(data = 0.,
                           nrow = dim(x)[1],
                           ncol = nobj)
        x <<- x
        for (i in 1:dim(x)[1]) {
          res <- func(i)
          newfeval[i, ] <- res[1:nobj]
          additional_eval <- rbind(additional_eval, res[(nobj + 1):length(res)])
        }
        nadditional <- ncol(additional_eval)
        if (length(res) < (nobj + 1)) {
          additional_eval <- NULL
          nadditional <- 0
        }
      } else {
        newfeval <- func(x)
        additional_eval <- NULL
        nadditional <- 0
      }

      # update the number of calls
      nrun <- nrun + dim(x)[1]

      # show progress
      if (verbose) setTxtProgressBar(pb, min(nrun, maxrun) / maxrun)

      # deal with NaN value
      detect_nan <- is.na(newfeval)
      set_ok <- !rowSums(detect_nan)
      if (length(set_ok[set_ok == TRUE]) == 0) {
        if (carallel == 1) stopCluster(cl)
        if (verbose) {
          message("Optimization failed")
          close(pb)
        }
        return(list("success" = FALSE, "message" = "No feasible points! Try to increase the size of the population..."))
      }
      newfeval <- newfeval[set_ok, ]
      dim(newfeval) <- c(sum(set_ok, na.rm = TRUE), nobj)
      x <- x[set_ok, ]
      dim(x) <- c(sum(set_ok, na.rm = TRUE), nvar)
      additional_eval <- additional_eval[set_ok, ]
      probj <- probj[set_ok, ]
      pop1 <- rbind(pop, cbind(x, newfeval, additional_eval))

      # decrease population
      matobj <- pop1[, (nvar + 1):(nvar + nobj)]
      dim(matobj) <- c(dim(pop1)[1], nobj)
      ind <- decrease_pop(matobj, minmax, prec, archsize, popsize)

      # archive
      arch <- matrix(pop1[ind$arch, ], nrow = length(ind$arch), ncol = nobj + nvar + nadditional)

      # population update
      pop <- pop1[c(ind$arch, ind$pop), ]
      dim(pop) <- c(length(c(ind$arch, ind$pop)), nvar + nobj)
      param_arch <- arch[, 1:nvar]
      crit_arch <- matrix(arch[, (nvar + 1):(nvar + nobj)], nrow = length(ind$arch), ncol = nobj)
      if (nadditional > 0) {
        additional_eval <- matrix(arch[, (nvar + nobj + 1):(nvar + nobj + nadditional)], nrow = length(ind$arch), ncol = nadditional)
      }

      # Records on criteria
      a <- c(lapply(c(1:nobj), function(i) max(crit_arch[, i])))
      maxcrit <- as.data.frame(a, col.names = objnames)
      a <- c(lapply(c(1:nobj), function(i) min(crit_arch[, i])))
      mincrit <- as.data.frame(a, col.names = objnames)
      crit <- mincrit
      crit[minmax] <- maxcrit[minmax]
      save_crit <- cbind(save_crit, c(nrun, t(crit)))

      # Graphs
      if (graph == TRUE) plot_population(MatObj = crit_arch, nobj, ngen, nrun, objnames, MatEvol = save_crit, popsize)

      # Online saves
      if (writefile == TRUE) {
        if (write_gen == TRUE) {
          listsave$pmt <- paste(listsave$RadPmt, "_gen", ngen, ".txt", sep = "")
          listsave$obj <- paste(listsave$RadObj, "_gen", ngen, ".txt", sep = "")
          if (ecrit_total_pop == 1) {
            listsave$totalpop <- paste(listsave$RadPop, "_gen", ngen, ".txt", sep = "")
          }
        }

        write.table(param_arch, listsave$pmt, row.names = FALSE, col.names = FALSE)
        write.table(cbind(crit_arch, additional_eval), listsave$obj, row.names = FALSE, col.names = FALSE)
        write.table(t(save_crit), listsave$evol, row.names = FALSE, col.names = FALSE)
        if (ecrit_total_pop == 1) {
          write.table(pop, listsave$totalpop, row.names = FALSE, col.names = FALSE)
        }
      }
    }

    # Close the progress bar
    if (verbose) close(pb)

    # Compute the first order derivatives
    if (sensitivity == TRUE) {
      if (verbose) message("Computing the sensitivity of the Pareto front...")
      dx <- 0.0001
      dxinv <- 1. / dx
      jacobian <- list()

      xopt <- param_arch
      nfront <- dim(crit_arch)[1]
      dim(xopt) <- c(nfront, nvar)
      x <- xopt
      for (k in 1:nobj) {
        nameId <- paste("Jacobian_", toString(k), sep = "")
        jacobian[[nameId]] <- matrix(data = 0., nrow = nfront, ncol = nvar)
        dim(jacobian[[nameId]]) <- c(nfront, nvar)  # even if only one variable it must be a matrix
      }
      for (j in 1:nvar) {
        x[, j] <- x[, j] + dx
        if (carallel == 1) {
          newfeval <- NULL
          clusterExport(cl = cl, varlist = c("x", "worklist"), envir = environment())
          res <- parLapply(cl, 1:nfront, func)
          for (e in 1:nfront) {
            newfeval <- rbind(newfeval, as.numeric(res[[e]][1:nobj]))
          }
        } else if (carallel == 0 | carallel == 2) { #sequential calls or user-defined choice
          newfeval <- matrix(data = 0.,
                             nrow = nfront,
                             ncol = nobj)
          x <<- x
          for (e in 1:nfront) {
            res <- func(e)
            newfeval[e, ] <- res[1:nobj]
          }
        }
        for (k in 1:nobj) {
          jacobian[[k]][, j] <- (newfeval[, k] - crit_arch[, k]) * dxinv
        }
        x[, j] <- xopt[, j]
      }
      nrun <- nrun + dim(x)[1] # update the number of calls
    }

    # Stop the // cluster
    if (carallel == 1) stopCluster(cl)

    end_time <- Sys.time()
    if (verbose) {
      message(paste("Done in", as.character(end_time - start_time), units(end_time - start_time), "-->", date()))
      message(paste("Size of the Pareto front :", as.character(dim(crit_arch)[1])))
      message(paste("Number of calls :", as.character(nrun)))
    }

    return(list(
      "success" = TRUE,
      "parameters" = param_arch,
      "objectives" = cbind(crit_arch, additional_eval),
      "derivatives" = jacobian,
      "save_crit" = t(save_crit),
      "total_pop" = pop,
      "gpp" = gpp
    ))
  }
