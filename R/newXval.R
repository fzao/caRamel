#' newXval
#' 
#' La fonction newXval genere une nouvelle population de jeux de parametres en suivant les cinq regles de caRamel.
#'  
#' @param param matrice [ Nvec , NPar ] des parametres de la population courante
#' @param crit matrice [ Nvec , NObj ] des criteres associes
#' @param isperf vecteur de booleens de longueur NObj, TRUE si maximisation de l'objectif, FALSE sinon
#' @param sp Variance a priori des parametres
#' @param bounds bornes inf et sup des parametres [ NPar , 2 ]
#' @param repart_gene : matrice de longueur 4 donnant le nombre de jeux a generer avec chaque regle : 1 Interpolation dans les simplexes du front, 2 Extrapolation selon les directions des aretes "orthogonales" au front, 3 Tirages aleatoires avec matrice de variance-covariance prescrite, 4 Recombinaison par blocs fonctionnels
#' @param blocs liste de vecteurs d'entier contenant les blocs fonctionnels de parametres
#' @param fireworks : booleen, TRUE si on teste une variation aleatoire sur chaque parametre et chaque maximum de F.O.
#' @return xnouv matrice de nouveaux vecteurs [ sum(Repart_Gene) + eventuellement (nobj+1)*nvar si fireworks , NPar ]
#' @return project_crit position supposee des nouveau vecteurs dans l'espace des criteres : [ sum(Repart_Gene)+ eventuellement (nobj+1)*nvar si fireworks , NObj ];
#' @author F. Zaoui
#' @export

newXval <-
  function(param,
           crit,
           isperf,
           sp,
           bounds,
           repart_gene,
           blocs,
           fireworks) {

    nobj <- dim(crit)[2]
    npar <- dim(param)[2]
    nvec <- dim(param)[1]
    a <- 3 / 8
    
    # Normalisation des differents objectifs
    obj <- matrix(0, nrow = dim(crit)[1], ncol = nobj)
    for (i in 1:nobj) {
      obj[, i] <- (val2rank(crit[, i], 3) - a) / (nvec + 1 - 2 * a)
    }
    obj[, !isperf] <- 1 - obj[, !isperf]
    
    Fo <- dominate(obj)
    param_arch <- as.matrix(param[Fo == 1, ])
    if (dim(param_arch)[2]<npar){param_arch <- t(param_arch)} # pour garder nb_param colonnes
    obj_arch <- obj[Fo == 1, ]
    crit_arch <- crit[Fo == 1, ]
    
    n_inter <- repart_gene[1] #Nombre maximum de nouveaux jeux crees par interpolation dans les simplexes
    n_extra <- repart_gene[2] #Nombre maximum de nouveaux jeux crees par extrapolation sur les directions des aretes
    n_cov <- repart_gene[3]   #Nombre maximum de nouveaux jeux crees par tirages aleatoires avec matrice de covariance prescrite
    n_recomb <- repart_gene[4]#Nombre maximum de nouveaux jeux crees par recombinaisons
    
    xnouv <- NULL
    project_crit <- NULL
    
    #**********************************************************************
    #            TRIANGULATION DE L'ESPACE DES OBJECTIFS (SI BESOIN)      *
    #**********************************************************************
    
    if (n_inter > 0 | n_extra > 0) {
      simplices <- delaunayn(obj)
      
      # Pour chaque simplexe on calcule le nombre nf de sommets appartenant au front de Pareto 
      nf <-
        apply((matrix(
          Fo[simplices],
          nrow = dim(simplices)[1],
          ncol = dim(simplices)[2]
        ) == 1), 1, sum)
      
      # On ne garde que les simplexes ayant au moins un sommet sur le front de Pareto
      ix <- which(nf > 0)
      simplices <- simplices[ix, ]
      simplices <- matrix(data = simplices,ncol=(nobj+1)) # Pour garder une matrice meme si un seul simplex
      nbsimp <- length(ix)
      
      # Decomposition en aretes, calcul du volume des simplexes conserves
      na <- nobj * (nobj + 1) / 2 #nombre d'aretes dans un simplexe en dimension NObj
      
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
    
    #**********************************************************************
    #                   GENERATION DES NOUVEAUX POINTS                    *
    #**********************************************************************
    
    # Nouveaux jeux crees par interpolation dans les simplexes  
    if (n_inter > 0) {
      carain <- Cinterpole(param, crit, simplices, volume, n_inter)
      xnouv <- rbind(xnouv, carain$xnouv)
      project_crit <- rbind(project_crit, carain$pcrit)
    }
    
    # Nouveaux jeux crees par extrapolation sur les directions des aretes
    if (n_extra > 0) {
      caraex <- Cextrapole(param, crit, oriedge, ledge, n_extra)
      xnouv <- rbind(xnouv, caraex$xnouv)
      project_crit <- rbind(project_crit, caraex$pcrit)
    }
    
    # Nouveaux jeux crees par tirages aleatoire avec matrice de covariance prescrite
    if (n_cov > 0) {
      # Calcul de la matrice de variance-covariance sur une sous-population de reference
      iref <- sort(unique(c(simplices))) # Population de reference : ensembles des points sommet d'un simplexe "frontal"
      xref <- as.matrix(param[iref, ])
      xcov <- Cusecovar(xref, sqrt(2), n_cov)
      critcov <- matrix(NaN, nrow = n_cov, ncol = nobj)
      xnouv <- rbind(xnouv, xcov)
      project_crit <- rbind(project_crit, critcov)
    }
    
    # Nouveaux jeux crees avec recombinaisons par blocs
    if (dim(param_arch)[1] > 1 & n_recomb > 0) {
      xrecomb <- Crecombination(param_arch, blocs, n_recomb)
      critrec <- matrix(NaN, nrow = n_recomb, ncol = nobj)
      xnouv <- rbind(xnouv, xrecomb)
      project_crit <- rbind(project_crit, critrec)
    }
    
    # Nouveaux jeux crees par variations independantes des parametres
    if (fireworks) {
      
      if (sum(Fo == 1) == 1){              # Dans le cas ou seulement 1 point sur le front
        obj_arch <- matrix(obj_arch,1,nobj)
        crit_arch <- matrix(crit_arch,1,nobj)
      }
      
      sp <- matrix(sp, nrow = 1, ncol = length(sp))
      
      #Points maximisant chaque F.O. individuellement
      m <- apply(obj_arch, 2, max)
      ipp <- apply(obj_arch, 2, which.max)
      
      # Maxi-min (point "central" du front)
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
