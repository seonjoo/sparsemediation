% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ginv.largep.R
\name{ginv.largep}
\alias{ginv.largep}
\title{Compute inverse, squareroot and inverse of the square root of the covariance}
\usage{
ginv.largep(x.c, sqrtmat = TRUE, sqrtinvmat = TRUE)
}
\arguments{
\item{x.c}{a data matrix (n x p)}

\item{sqrtmat}{(default = TRUE) return square root of the covariance matrix}

\item{sqrtinvmat}{(default = TRUE) return invere of the squared root of the covariance matrix}
}
\value{
squareroot of matrix
}
\description{
This is a miscelaneous function.
}
\examples{
set.seed(1234)
x=scale(matrix(100*20,100,20),center=TRUE, scale=FALSE)
ginv.largep(x)

}
\references{
TBA
}
\author{
Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
}
\keyword{L1penalization}
\keyword{highdimensional}
\keyword{mediation}
