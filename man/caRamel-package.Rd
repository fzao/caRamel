\name{caRamel-package}
\alias{caRamel-package}
\docType{package}

\title{caRamel optimizer}

\description{Automatic Calibration by Evolutionary Multi Objective Algorithm}

\details{
\emph{caRamel} is a package for multi-objective optimization of complex environmental models.

The algorithm is a hybrid of the \emph{MEAS} algorithm (Efstratiadis and Koutsoyiannis, 2005)
by using the directional search method based on the simplexes of the objective space 
and the \emph{epsilon-NGSA-II} algorithm with the method of classification of the parameter vectors archiving management by epsilon-dominance
(Reed and Devireddy, 2004).

The main function of the package is \emph{caRamel()}.

This function uses all the other functions of the package.

An example of an hydrological optimization is available on the following presentation: \href{https://www.researchgate.net/publication/341931417_Automatic_Calibration_by_Evolutionary_Multi_Objective_Algorithm_the_caRamel_R_package}{useR! 2019}

}
\author{Fabrice Zaoui, Nicolas Le Moine, Celine Monteil (EDF R&D - LNHE)}
\references{
   
   Efstratiadis, A. and Koutsoyiannis, D. (2005) \emph{The multi-objective evolutionary annealing-simplex method and its application in calibration hydrological models},
   in EGU General Assembly 2005, Geophysical Research Abstracts, Vol. 7, Vienna, 04593, European Geophysical Union. doi:10.13140/ RG.2.2.32963.81446.
   
   Le Moine, N. (2009) \emph{Description d’un algorithme génétique multi-objectif pour la calibration d’un modèle pluie-débit} (in French).
   Post-Doctoral Status Rep. 2, UPMC/EDF, 13 pp.
   
   Reed, P. and Devireddy, D. (2004) \emph{Groundwater monitoring design: a case study combining epsilon-dominance archiving and automatic parameterization for the NSGA-II},
   in Coello-Coello C, editor. Applications of multi-objective evolutionary algorithms, Advances in natural computation series,
   vol. 1, pp. 79-100, World Scientific, New York. doi:10.1142/9789812567796_0004.
}
\keyword{package}

