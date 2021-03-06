#' @title Sparse mediation for high-dimensional mediators
#'
#' Fit a mediation model via penalized maximum likelihood and structural equation model.
#' The regularization path is computed for the lasso or elasticnet penalty at a grid of
#' values for the regularization parameter lambda. Currently, mediation analysis is developed based on gaussian assumption.
#'
#' Multiple Mediaton Model:
#' (1) M = Xa + e1
#' (2) Y = Xc' + Mb + e2
#' And in the optimization, we do not regularize c', due to the assumption of partial mediation.
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param tol (default -10^(-5)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:50)/125)) tuning parameter for L1 penalization
#' @param lambda2 (default=c(0.2,0.5)) tuning parameter for L1 penalization for covariance matrix, used only for p>n.
#' @param alpha (default=1) tuning parameter for L2 penalization
#' @param tau (default=1) tuning parameter for differentail weight between paths a (X -> M) and b (M -> Y)
#' @param verbose (default=TRUE) print progress.
#' @param Omega.out (defult=TRUE) output Omega estimates
#' @return c: directeffect per each tuning parameter lambda. length(lambda)-dimensional vector
#' @return hatb: Path b (M->Y given X) estimates: V-by-lenbth(lambda)  matrix
#' @return hata: Path a (X->M) estimates: V-by-lenbth(lambda)  matrix
#' @return medest: Mediation estimates (a*b): V-by-lenbth(lambda)  matrix
#' @return alpha: a scalor of the numing parameter for L2 regularization
#' @return lambda: a vector of tuning parameters for L1-penalization
#' @return tau: weight used.
#' @return nump: Number of selected mediation paths
#' @return Omega Estimated covariance matrix of the mediator
#' @examples
#' library(sparsemediation)
#' N=100
#' V=50
#' set.seed(1234)
#' a = rep(0,V);a[1:3]<-5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' sparse.mediation(X,M,Y)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import glmnet
#' @import Matrix
#' @import QUIC
#' @importFrom stats var predict
#' @export

sparse.mediation = function(X,M,Y,tol=10^(-5),max.iter=50,
                            lambda = log(1+(1:30)/100),
                            lambda2 = c(0.3,0.5),
                            alpha=1,tau=1,verbose=FALSE,Omega.out=FALSE){

  re=new.env()
  re=as.list(re)
  V = ncol(M)
  N=nrow(M)

 # if(N > 2*V){
#    re=sparse.mediation.old(X,M,Y,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha,tau=tau,verbose=verbose)
 # }else{
  re=sparse.mediation.largep_omega(X,M,Y,tol=tol,max.iter=max.iter,lambda1 = lambda,
                                                lambda2=lambda2,alpha=alpha,tau=tau,verbose=verbose,
                                   Omega.out=Omega.out)
  #}
  return(re)
}
