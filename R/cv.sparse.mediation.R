#' Conduct K-fold cross validation for sparse mediation with elastic net with multiple tuning parameters
#'
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param K (default=5) number of cross-validation folds
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:30)/100)) tuning parameter for L1 penalization
#' @param lambda2 (default=c(0.2,0.5)) tuning parameter for inverse covariance matrix sparsity. Only used if n>(2*V).
#' @param alpha (defult=1) tuning parameter for L2 penalization
#' @param tau (default=1) tuning parameter for differential weight for L1 penalty.
#' @param glmnet.penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @param multicore (default=1) number of multicore
#' @param seednum (default=10000) seed number for cross validation
#' @return cv.lambda: optimal lambda
#' @return cv.tau: optimal tau
#' @return cv.alpha: optimal tau
#' @return cv.mse: minimum MSE value
#' @return mse: Array of MSE, length(alpha) x length(lambda) x length (tau)
#' @return lambda: vector of lambda
#' @return tau: vector of tau used
#' @return alpha: vector of alpha used
#' @return z: cross-valication results
#' @examples
#' N=200
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  as.vector(X + M %*% b + rnorm(N))
#' system.time(cvfit<-cv.sparse.mediation(X, M, Y, tol = 10^(-10), K = 4, max.iter = 100,lambda = log(1 + (1:10)/25), alpha = 1, tau=c(0.25,0.5,1,2,4,8),multicore = 4, seednum = 1e+06))
#' cvfit$cv.lambda
#' cvfit$cv.tau
#' fit<-sparse.mediation(X,M,Y,tol=10^(-10),max.iter=100,lambda = cvfit$cv.lambda, alpha=cvfit$cv.alpha,tau=cvfit$cv.tau)
#' nonzerolist = c(0,as.numeric(abs(c(fit[[1]]$hatb, fit[[1]]$hata))<0.0001))
#' refit=sparse.mediation.old(X,M,Y,lambda = 100000,alpha=cvfit$cv.alpha, glmnet.penalty.factor=nonzerolist)
#' cbind(a*b,refit$medest)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca glmnet
#' @import parallel
#' @import MASS
#' @import glmnet
#' @export


cv.sparse.mediation= function(X,M,Y,tol=10^(-10),K=5,max.iter=100,
                              lambda = log(1+(1:15)/50),
                              lambda2 = c(0.2,0.5),
                              alpha=1,tau=c(0.5,1,2),
                              multicore=1,seednum=1000000,verbose=FALSE){

  ## Center all values
  N = nrow(M)
  V = ncol(M)
  Y.mean=mean(Y)
  X.mean=mean(X)
  M.mean=apply(M,2,mean)
  Y.sd=sqrt(var(Y))
  X.sd=sqrt(var(X))
  M.sd=sqrt(apply(M,2,var))


  Y = scale(Y,center=TRUE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=TRUE,scale=TRUE)

  ###K-fold cross-validation
  set.seed(seednum)
  cvid = (rep(1:K, each=ceiling(N/K))[1:N])[sort.int(rnorm(N),index.return=TRUE)$ix]

  if(multicore>1){
    z<-mclapply(1:K, function(fold){
      re=sparse.mediation(Y=Y[cvid!=fold,],X= X[cvid!=fold,],M=M[cvid!=fold,], lambda=lambda, lambda2=lambda2,tol=tol,alpha=alpha,
                       tau=tau,verbose=verbose)
      y=Y[cvid==fold,];x=X[cvid==fold,] ;m=M[cvid==fold,]
      mse=unlist(lapply(1:ncol(re$hata), function(loc){ mean((y - re$c[loc] - m %*%re$hatb[,loc])^2) + mean((m - x %*% t(re$hata[,loc]))^2)}))
      return(list(re=re,mse=mse))},mc.cores=multicore)
  }else{
    z<-lapply(1:K, function(fold){
      re=sparse.mediation(Y=Y[cvid!=fold,],X= X[cvid!=fold,],M=M[cvid!=fold,], lambda=lambda, lambda2=lambda2,tol=tol,alpha=alpha,
                          tau=tau,verbose=verbose)
      y=Y[cvid==fold,];x=X[cvid==fold,] ;m=M[cvid==fold,]
      mse=unlist(lapply(1:ncol(re$hata), function(loc){ mean((y - re$c[loc] - m %*%re$hatb[,loc])^2) + mean((m - x %*% t(re$hata[,loc]))^2)}))
      return(list(re=re,mse=mse))})
  }
  mseest=apply(do.call(cbind,lapply(z,function(x)x$mse)),1,mean)
  minloc=which.min(mseest)
  min.lambda1=z[[1]]$re$lambda1[minloc]
  min.lambda2=z[[1]]$re$lambda2[minloc]
  min.alpha=z[[1]]$re$alpha[minloc]
  min.tau=z[[1]]$re$tau[minloc]

  return(list(cv.lambda1=min.lambda1,cv.lambda2=min.lambda2,cv.tau=min.tau, cv.alpha=min.alpha,cv.mse=mseest[minloc],mse=mseest))
}


