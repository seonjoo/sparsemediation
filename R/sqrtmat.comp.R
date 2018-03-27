#' Compute squareroot matrix
#' @param mat One-dimensional predictor
#' @param threshold
#' @return squareroot of matrix
#' @examples
#' set.seed(1234)
#' x=matrix(100*20,100,20)
#' covmat=cov(x)
#' sqrtmat.comp(covmat)
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation L1penalization
#' @import parallel
#' @import MASS
#' @import glmnet
#' @export

sqrtmat.comp<-function(mat,thresh=10^(-20)){
  eigenmat=svd(mat)
  #  print(paste('Dimension of mat:', dim(mat)))
  ncomp=sum(eigenmat$d>thresh)
  #print(ncomp)
  # print(eigenmat$d)
  if (ncomp<2){
    sqmat=as.matrix(eigenmat$v[,1]) %*% sqrt(eigenmat$d[1]) %*% t(as.matrix(eigenmat$v[,1]))
  }else{sqmat=eigenmat$v[,1:ncomp] %*% diag(sqrt(eigenmat$d[1:ncomp])) %*% t(eigenmat$v[,1:ncomp])
  }


  return(sqmat)
}

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

ginv.largep0<-function(x,sqrtmat=TRUE, sqrtinvmat=TRUE){
  n=nrow(x)
  x.c=scale(x, center=TRUE, scale=FALSE)
  xxt.inv= ginv( x.c %*% t(x.c)/(n-1))
  tmp = xxt.inv %*% x.c/(sqrt(n-1))
  sqrt.mat=sqrt.invmat=NULL
  if (sqrtinvmat==TRUE){
    sqrt.mat=t(sqrtmat.comp(xxt.inv) %*% x.c/(sqrt(n-1))) %*% x.c/sqrt(n-1)
  }
  if (sqrtinvmat==TRUE){
    sqrt.invmat=t(sqrtmat.comp(xxt.inv) %*% x.c/(sqrt(n-1))) %*% xxt.inv %*% x.c/sqrt(n-1)
  }
  return(list(inv=t(tmp) %*% tmp, sqrtinv=sqrt.invmat, sqrtmat=sqrt.mat))
}
