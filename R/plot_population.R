#' Plotting of a population of objectives
#'
#' Plot graphs the population regarding each couple of objectives
#'
#' @param MatObj : matrix of the objectives [NInd, nobj]
#' @param nobj : number of objectives
#' @param ngen : number of generations (optional)
#' @param nrun : number of model evaluations (optional)
#' @param objnames : vector of objectives names (optional)
#' @param MatEvol : matrix of the evolution of the optimal objectives (optional)
#' @param popsize : integer, size of the initial population (optional)
#'
#' @examples
#' # Definition of the population
#' Pop <- matrix(runif(300), 100, 3)
#' # Call the function
#' plot_population(MatObj = Pop, nobj = 3, objnames = c("Obj1", "Obj2", "Obj3"))
#'
#' @author Celine Monteil

plot_population <- function(MatObj, nobj, ngen = NULL, nrun = NULL, objnames = NULL, MatEvol = NULL, popsize = 0) {

  # Graphs
  info <- paste("ngen=", ngen, ", nrun=", nrun, sep = "")
  nb_fen <- choose(n = nobj, k = 2) + 1
  if (nb_fen <= 4) {
    par(mfrow = c((nobj - 1), 2))
  } else {
    par(mfrow = c(3, floor(nb_fen / 3) + 1))
  }
  l <- seq(1, nobj)
  for (i_fig in 1:(nobj - 1)) {
    l_tmp <- l[-i_fig]
    l_tmp <- l_tmp[l_tmp > i_fig]
    for (i_fig2 in l_tmp) {
      plot(MatObj[, i_fig], MatObj[, i_fig2], xlab = objnames[i_fig], ylab = objnames[i_fig2])
    }
  }

  if (!is.null(MatEvol)) {
    plot(x = MatEvol[1, ], y = MatEvol[2, ],
         ylim = c(min(MatEvol[-1, ]), max(MatEvol[-1, ])),
         xlab = info, ylab = "Optimal Criteria", pch = 0)

    lapply(c(2:nobj), function(i) {
       points(x = MatEvol[1, ], y = MatEvol[i + 1, ], col = i, pch = (i - 1))
    })

    xlgd <- popsize + (nrun - popsize) * 2 / 3
    ylgd <- min(MatEvol[-1, ]) + (max(MatEvol[-1, ]) - min(MatEvol[-1, ])) / 2
    legend(xlgd, ylgd, legend = objnames, col = 1:nobj, pch = 0:(nobj - 1))
  }
}
