% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dominated.R
\name{dominated}
\alias{dominated}
\title{Rows domination of a matrix by a vector}
\usage{
dominated(x, Y)
}
\arguments{
\item{x}{: row vecteur}

\item{Y}{: matrix}
}
\value{
D : vector of booleans
}
\description{
indicates which rows of the matrix Y are dominated by the vector (row) x
}
\examples{
# Definition of the parameters
Y <- matrix(rexp(200), 100, 2)
x <- Y[1,]
# Call the function
res <- dominated(x, Y)

}
\author{
Alban de Lavenne, Fabrice Zaoui
}
