#' Conduct sparse mediation with elastic net with multiple tuning parameters
#'
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:30)/100)) tuning parameter for L1 penalization
#' @param alpha (defult=c(1:4)/4) tuning parameter for L2 penalization
#' @param figure (defult=NULL) print figures for mean predictive errors by tuning parameters alpha and lambda
#' @param glmnet.penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @return re list of sparse.mediation per each alpha
#' @examples
#' N=100
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' sparse.mediation.twoway(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:25)/50))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca glmnet
#' @export

sparse.mediation.twoway = function(X,M,Y,tol=10^(-10),max.iter=100,
                                   lambda = log(1+(1:30)/100),alpha=(1:4)/4,
                                   glmnet.penalty.factor=c(0,rep(1,2*V))){
  library(parallel)
  library(MASS)
  library(glmnet)
  re=new.env()
  re=as.list(re)
  for (j in 1:length(alpha)){
    re[[j]]=sparse.mediation(X,M,Y,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha[j],glmnet.penalty.factor=glmnet.penalty.factor)
    re[[j]]$alpha=alpha[j]
  }
  return(re)
}

