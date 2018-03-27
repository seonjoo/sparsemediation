#' Conduct sparse mediation with elastic net
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
#' @param alpha (default=1) tuning parameter for L2 penalization
#' @param tau (default=1) tuning parameter for differentail weight between paths a (X -> M) and b (M -> Y)
#' @param figure (defult=NULL) print figures for mean predictive errors by tuning parameters alpha and lambda
#' @param glmnet.penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @return c: directeffect per each tuning parameter lambda. length(lambda)-dimensional vector
#' @return hatb: Path b (M->Y given X) estimates: V-by-lenbth(lambda)  matrix
#' @return hata: Path a (X->M) estimates: V-by-lenbth(lambda)  matrix
#' @return medest: Mediation estimates (a*b): V-by-lenbth(lambda)  matrix
#' @return alpha: a scalor of the numing parameter for L2 regularization
#' @return lambda: a vector of tuning parameters for L1-penalization
#' @return tau: weight used.
#' @return nump: Number of selected mediation paths
#' @examples
#' N=100
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' sparse.mediation(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:25)/50))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import glmnet
#' @export

sparse.mediation = function(X,M,Y,tol=10^(-5),max.iter=50,
                                   lambda = log(1+(1:30)/100),alpha=1,tau=1){
#  library(parallel)
#  library(MASS)
#  library(glmnet)


  re=new.env()
  re=as.list(re)
  V = ncol(M)
  weights = c(0,rep(1,V), tau*rep(1,V))
  for (j in 1:length(alpha)){
    if(N > V){
    re[[j]]=sparse.mediation.old(X,M,Y,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha[j],glmnet.penalty.factor=weights)
    }else{re[[j]]=sparse.mediation.largep(X,M,Y,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha[j],glmnet.penalty.factor=weights)}

    re[[j]]$alpha=alpha[j]
    re[[j]]$tau=tau
  }
  return(re)
}
