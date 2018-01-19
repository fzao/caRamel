#' newXval
#' 
#' generates a new population of parameter sets following the five rules of caRamel
#'  
#' @param param : matrix [ Nvec , NPar ] of parameters of the current population
#' @param crit : matrix [ Nvec , NObj ] of associated criteria
#' @param isperf : vector of Booleans of length NObj, TRUE if maximization of the objective, FALSE otherwise
#' @param sp : variance a priori of the parameters
#' @param bounds : lower and upper bounds of parameters [ NPar , 2 ]
#' @param repart_gene : matrix of length 4 giving the number of games to be generated with each rule: 1 Interpolation in the simplexes of the front, 2 Extrapolation according to the directions of the edges "orthogonal" to the front, 3 Random draws with prescribed variance-covariance matrix, 4 Recombination by functional blocks
#' @param blocks : list of integer vectors containing function blocks of parameters
#' @param fireworks : boolean, TRUE if one tests a random variation on each parameter and each maximum of O.F.
#' @return xnouv : matrix of new vectors [ sum(Repart_Gene) + eventually (nobj+1)*nvar if fireworks , NPar ]
#' @return project_crit: assumed position of the new vectors in the criteria space: [ sum(Repart_Gene)+ eventually (nobj+1)*nvar if fireworks , NObj ];
#' @author Fabrice Zaoui

newXval <-
  function(param,
           crit,
           isperf,
           sp,
           bounds,
           repart_gene,
           blocks,
           fireworks) {

    nobj <- dim(crit)[2]
    npar <- dim(param)[2]
    nvec <- dim(param)[1]
    a <- 3 / 8
    
    # Standardization of the different objectives
    obj <- matrix(0, nrow = dim(crit)[1], ncol = nobj)
    for (i in 1:nobj) {
      obj[, i] <- (val2rank(crit[, i], 3) - a) / (nvec + 1 - 2 * a)
    }
    obj[, !isperf] <- 1 - obj[, !isperf]
    
    Fo <- dominate(obj)
    param_arch <- as.matrix(param[Fo == 1, ])
    if (dim(param_arch)[2]<npar){param_arch <- t(param_arch)} # to keep nb_param columns
    obj_arch <- obj[Fo == 1, ]
    crit_arch <- crit[Fo == 1, ]
    
    n_inter <- repart_gene[1]  # Maximum number of new sets created by interpolation in simplexes
    n_extra <- repart_gene[2]  # Maximum number of new sets created by extrapolation on the directions of the edges
    n_cov <- repart_gene[3]    # Maximum number of new sets created by random draws with prescribed covariance matrix
    n_recomb <- repart_gene[4] # Maximum number of new games created by recombinations
    
    xnouv <- NULL
    project_crit <- NULL
    
    #**************************************************************
    #       TRIANGULATION OF SPACE OF OBJECTIVES (IF NEEDED)      *
    #**************************************************************
    
    if (n_inter > 0 | n_extra > 0) {
      simplices <- delaunayn(obj)
      
      # For each simplex calculation of the number of vertices belonging to the Pareto front
      nf <-
        apply((matrix(
          Fo[simplices],
          nrow = dim(simplices)[1],
          ncol = dim(simplices)[2]
        ) == 1), 1, sum)
      
      # Keep only the simplexes with at least one vertex on the Pareto front
      ix <- which(nf > 0)
      simplices <- simplices[ix, ]
      simplices <- matrix(data = simplices,ncol=(nobj+1)) # To keep a matrix even if a single simplex
      nbsimp <- length(ix)
      
      # Decomposition in edges, calculation of the volume of the simplexes preserved
      na <- nobj * (nobj + 1) / 2 #number of edges in a NObj size simplex
      
      volume <- matrix(0, nrow = nbsimp, ncol = 1)
      oriedge <- NULL
      ledge <- NULL
      
      for (s in 1:nbsimp) {
        P <- simplices[s, ]
        S <- obj[P, ]
        volume[s] <- vol_splx(S)
        direc <- Dimprove(S, Fo[P])
        if (!is.null(direc$oriedge)) {
          oriedge <-
            rbind(oriedge,  matrix(P[direc$oriedge], nrow = dim(direc$oriedge)[1]))
          ledge <- c(ledge, direc$ledge)
        }
      }
      
      unik <- !duplicated(oriedge)
      #ix <- seq_along(oriedge)[unik]
      # tmp <- oriedge[unik]
      oriedge <- oriedge[unik, ]
      ledge <- ledge[unik]
    }
    
    #**************************************************************
    #                 GENERATION OF NEW POINTS                    *
    #**************************************************************
    
    # New sets created by interpolation in simplexes  
    if (n_inter > 0) {
      carain <- Cinterp(param, crit, simplices, volume, n_inter)
      xnouv <- rbind(xnouv, carain$xnouv)
      project_crit <- rbind(project_crit, carain$pcrit)
    }
    
    # New sets created by extrapolation on the directions of the edges
    if (n_extra > 0) {
      caraex <- Cextrap(param, crit, oriedge, ledge, n_extra)
      xnouv <- rbind(xnouv, caraex$xnouv)
      project_crit <- rbind(project_crit, caraex$pcrit)
    }
    
    # New sets created by random draws with prescribed covariance matrix
    if (n_cov > 0) {
      # Calculation of the variance-covariance matrix on a reference subpopulation
      iref <- sort(unique(c(simplices))) # Reference population: all vertex points of a "frontal" simplex
      xref <- as.matrix(param[iref, ])
      xcov <- Cusecovar(xref, sqrt(2), n_cov)
      critcov <- matrix(NaN, nrow = n_cov, ncol = nobj)
      xnouv <- rbind(xnouv, xcov)
      project_crit <- rbind(project_crit, critcov)
    }
    
    # New sets created with block recombinations
    if (dim(param_arch)[1] > 1 & n_recomb > 0) {
      xrecomb <- Crecombination(param_arch, blocks, n_recomb)
      critrec <- matrix(NaN, nrow = n_recomb, ncol = nobj)
      xnouv <- rbind(xnouv, xrecomb)
      project_crit <- rbind(project_crit, critrec)
    }
    
    # New sets created by independent variations of parameters
    if (fireworks) {
      if (sum(Fo == 1) == 1){              # In the case or only 1 point on the front
        obj_arch <- matrix(obj_arch,1,nobj)
        crit_arch <- matrix(crit_arch,1,nobj)
      }
      
      sp <- matrix(sp, nrow = 1, ncol = length(sp))
      
      #Points maximizing each O.F. individually
      m <- apply(obj_arch, 2, max)
      ipp <- apply(obj_arch, 2, which.max)
      
      # Maxi-min ("central" point of the front)
      m <- max(apply(obj_arch, 1, min))
      maximin <- which.max(apply(obj_arch, 1, min))
      ipp <- c(ipp, maximin)
      
      for (i in 1:length(ipp)) {
        xcloud <-
          matrix(param_arch[ipp[i], ],
                 nrow = npar,
                 ncol = npar,
                 byrow = TRUE)
        ccloud <-
          matrix(crit_arch[ipp[i], ],
                 nrow = npar,
                 ncol = nobj,
                 byrow = TRUE)
        devi <- c(rnorm(npar) * sp)
        xcloud <-
          xcloud + diag(x = devi,
                        nrow = length(devi),
                        ncol = length(devi))
        xcloud <- xcloud[sp > 0, ]
        ccloud <- ccloud[sp > 0, ]
        xnouv <- rbind(xnouv, xcloud)
        project_crit <- rbind(project_crit, ccloud)
      }
      
    }
    
    minp <-
      matrix(bounds[, 1],
             nrow = dim(xnouv)[1],
             ncol = npar,
             byrow = TRUE)
    maxp <-
      matrix(bounds[, 2],
             nrow = dim(xnouv)[1],
             ncol = npar,
             byrow = TRUE)
    ll = which(xnouv > maxp, arr.ind = TRUE)
    nelem <- dim(ll)[1]
    if (nelem > 0) {
      for (i in 1:nelem) {
        rnum <- ll[i, 1]
        cnum <- ll[i, 2]
        xnouv[rnum, cnum] <- maxp[rnum, cnum]
      }
    }
    ll = which(xnouv < minp, arr.ind = TRUE)
    nelem <- dim(ll)[1]
    if (nelem > 0) {
      for (i in 1:nelem) {
        rnum <- ll[i, 1]
        cnum <- ll[i, 2]
        xnouv[rnum, cnum] <- minp[rnum, cnum]
      }
    }
    
    return(list("xnouv" = xnouv, "pcrit" = project_crit))
  }
