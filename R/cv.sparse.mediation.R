#' Conduct K-fold cross validation for sparse medication with elastic net with multiple tuning parameters
#'
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param K (default=5) number of cross-validation folds
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:30)/100)) tuning parameter for L1 penalization
#' @param alpha (defult=c(1:4)/4) tuning parameter for L2 penalization
#' @param figure (defult=NULL) print figures for mean predictive errors by tuning parameters alpha and lambda
#' @param glmnet.penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @param multicore (default=1) number of multicore
#' @param seednum (default=10000) seed number for cross validation
#' @param disply devault=FALSE
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
#' system.time(cvfit<-cv.sparse.mediation(X, M, Y, tol = 10^(-10), K = 8, max.iter = 100,lambda = log(1 + (1:10)/25), alpha = 1, tau=c(0.25,0.5,1,2,4,8),figure = NULL, multicore = 8, seednum = 1e+06, display = TRUE))
#' cvfit$cv.lambda
#' cvfit$cv.tau
#' fit<-sparse.mediation(X,M,Y,tol=10^(-10),max.iter=100,lambda = cvfit$cv.lambda, alpha=cvfit$cv.alpha,tau=cvfit$cv.tau)
#' nonzerolist = c(0,as.numeric(abs(c(fit[[1]]$hatb, fit[[1]]$hata))<0.0001))
#' refit=sparse.mediation(X,M,Y,lambda = 100000,alpha=cvfit$cv.alpha, glmnet.penalty.factor=nonzerolist)
#' cbind(a*b,refit$medest)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca glmnet
#' @import parallel
#' @import MASS
#' @import glmnet
#' @export


cv.sparse.mediation= function(X,M,Y,tol=10^(-10),K=5,max.iter=100,
                              lambda = log(1+(1:15)/50),alpha=(0:4)/4,tau=c(0.5,1,2),
                              multicore=1,seednum=1000000){
  library(parallel)
  library(MASS)
  library(glmnet)
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
    options(cores = multicore)
    z<-mclapply(1:K, function(fold){cv.sparse.mediation.threeway.subroutine(fold, Y,X,M,cvid,lambda, alpha,tau,max.iter, tol)}, mc.cores=multicore)
  }else{
    z<-lapply(1:K, function(fold){cv.sparse.mediation.threeway.subroutine(fold, Y,X,M,cvid,lambda,alpha,tau,max.iter, tol)})
  }

  mseest=alphaest=lambdaest=tauest=array(NA,dim=c(length(alpha),length(lambda),length(tau)))


  for (k in 1:length(tau)){
    for (l in 1:length(alpha)){
    tmpmse=matrix(NA,K,length(lambda))
      for (j in 1:K){
        tmpmse[j,]=(z[[j]]$mse)[[k]][[l]]$mse
      }
      mseest[l,,k]=apply(tmpmse,2,sum)
      alphaest[l,,k]=rep(z[[j]]$mse[[k]][[l]]$alpha,length(lambda))
      lambdaest[l,,k]=z[[j]]$mse[[k]][[l]]$lambda
      tauest[l,,k]<-tau[k]

    }
  }
  minloc=which.min(mseest)
  min.lambda=lambdaest[minloc]
  min.alpha=alphaest[minloc]
  min.tau=tauest[minloc]

  return(list(cv.lambda=min.lambda,cv.tau=min.tau, cv.alpha=min.alpha,cv.mse=mseest[minloc],mse=mseest, lambda=lambdaest, tau=tauest,alpha=alphaest,z=z))

}

cv.sparse.mediation.twoway.subroutine<- function(k, Y,X,M,cvid,lambda,alpha,max.iter=100, tol=10^(-10)){
  Y.test=Y[cvid==k]#as.vector(scale(Y[cvid==k],center=TRUE,scale=FALSE))
  Y.train=Y[cvid!=k]
  X.test=matrix(X[cvid==k,],sum(cvid==k),1)#matrix(scale(X[cvid==k,],center=TRUE,scale=FALSE),sum(cvid==k),1)
  X.train=X[cvid!=k,]
  M.test=M[cvid==k,]#scale(M[cvid==k,],center=TRUE,scale=FALSE)
  M.train=M[cvid!=k,]
  #boxplot(t(alpha.train))
  r.train=sparse.mediation.twoway(X.train,M.train,Y.train,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha)
  mse = lapply(r.train, function(obj){cv.sparse.mediation.msecomputing(obj, Y.test, X.test, M.test)})
  return(list(r.train=r.train,mse=mse))
}

cv.sparse.mediation.threeway.subroutine<- function(fold, Y,X,M,cvid,lambda,alpha,tau,max.iter=100, tol=10^(-10)){
  Y.test=Y[cvid==fold]#as.vector(scale(Y[cvid==k],center=TRUE,scale=FALSE))
  Y.train=Y[cvid!=fold]
  X.test=matrix(X[cvid==fold,],sum(cvid==fold),1)#matrix(scale(X[cvid==k,],center=TRUE,scale=FALSE),sum(cvid==k),1)
  X.train=X[cvid!=fold]
  M.test=M[cvid==fold,]#scale(M[cvid==k,],center=TRUE,scale=FALSE)
  M.train=M[cvid!=fold,]
  #boxplot(t(alpha.train))
  r.train=sparse.mediation.threeway(X.train,M.train,Y.train,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha,tau=tau)
#  return(r.train)
#}
  mse<-c()
  length.alpha=length(alpha)
  for (l in 1:length(tau)){
    mse[[l]] = lapply(r.train[[(l-1)*length.alpha + (1:length.alpha)]], function(obj){cv.sparse.mediation.msecomputing(obj, Y.test, X.test, M.test)})
  }
  return(list(r.train=r.train,mse=mse))
}

cv.sparse.mediation.msecomputing<-function(obj, Y.test, X.test, M.test){
  V = (nrow(obj$beta))
  a.train=obj$hata
  b.train=obj$hatb
  c.train=matrix(obj$c,nrow=1)
  #	print(obj)
  yhat = X.test %*% c.train + M.test %*% b.train
  mse.m = rep(0,length(obj$lambda))
  for (j in 1:length(obj$lambda)){
    mhat=X.test %*% t(a.train[,j])
    mse.m[j]=sum((M.test - mhat)^2)
  }
  mse=(apply(Y.test-yhat, 2,function(x){sum(x^2)})) + mse.m
  if (length(obj$alpha)>0){
    return(list(mse=mse,lambda=obj$lambda,alpha=obj$alpha))
  }else{
    return(list(mse=mse,lambda=obj$lambda,alpha=NULL))
  }
}
