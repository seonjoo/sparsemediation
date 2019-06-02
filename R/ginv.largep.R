#' Compute inverse, squareroot and inverse of the square root of the covariance
#'
#' This is a miscelaneous function.
#'
#' @param x.c a data matrix (n x p)
#' @param sqrtmat (default = TRUE) return square root of the covariance matrix
#' @param sqrtinvmat (default = TRUE) return invere of the squared root of the covariance matrix
#' @return squareroot of matrix
#' @examples
#' set.seed(1234)
#' x=scale(matrix(100*20,100,20),center=TRUE, scale=FALSE)
#' ginv.largep(x)
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation L1penalization
#' @import parallel
#' @import MASS
#' @importFrom stats var predict
#' @export

ginv.largep<-function(x.c,sqrtmat=TRUE, sqrtinvmat=TRUE){
  #n=nrow(x)
  #x.c=scale(x, center=TRUE, scale=FALSE)
  xxt.inv= ginv( x.c %*% t(x.c))
  tmp = xxt.inv %*% x.c
  sqrt.mat=sqrt.invmat=NULL
  if (sqrtinvmat==TRUE){
    sqrt.mat=t(sqrtmat.comp(xxt.inv) %*% x.c) %*% x.c
  }
  if (sqrtinvmat==TRUE){
    sqrt.invmat=t(sqrtmat.comp(xxt.inv) %*% x.c) %*% xxt.inv %*% x.c
  }
  return(list(inv=t(tmp) %*% tmp, sqrtinv=sqrt.invmat, sqrtmat=sqrt.mat))
}
