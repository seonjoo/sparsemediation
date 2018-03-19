#' Conduct K-fold cross validation for sparse mediation with elastic net with multiple tuning parameters (old version)
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
#' @return re list of sparse.mediation per each alpha
#' @examples
#' N=100
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' cv.sparse.mediation.old(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:15)/50))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca glmnet



cv.sparse.mediation.old= function(X,M,Y,tol=10^(-10),K=5,max.iter=100,
                              lambda = log(1+(1:15)/50),alpha=(0:4)/4,
                              figure=NULL,multicore=1,seednum=1000000,display=TRUE){
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
    z<-mclapply(1:K, function(k){cv.sparse.mediation.twoway.subroutine(k, Y,X,M,cvid,lambda, alpha,max.iter, tol)}, mc.cores=multicore)
  }else{
    z<-lapply(1:K, function(k){cv.sparse.mediation.twoway.subroutine(k, Y,X,M,cvid,lambda,alpha,max.iter, tol)})
  }

  mseest=alphaest=lambdaest=matrix(NA,length(alpha),length(lambda))
  for (l in 1:length(alpha)){
    tmpmse=matrix(NA,K,length(lambda))
    for (j in 1:K){
      tmpmse[j,]=(z[[j]]$mse)[[l]]$mse
    }
    mseest[l,]=apply(tmpmse,2,sum)
    alphaest[l,]=rep((z[[j]]$mse)[[l]]$alpha,length(lambda))
    lambdaest[l,]=(z[[j]]$mse)[[l]]$lambda

  }
  minloc=which.min(mseest)
  min.lambda=lambdaest[minloc]
  min.alpha=alphaest[minloc]

  if (display==TRUE){
    colnum=rainbow(length(alpha),end=0.7)
    plot(lambdaest[1,],mseest[1,],type='l',col=colnum[1],ylim=range(mseest));points(lambdaest[1,],mseest[1,],pch=15,col=colnum[1])
    for (j in 1:length(alpha)){lines(lambdaest[j,],mseest[j,],type='l',col=colnum[j]);points(lambdaest[j,],mseest[j,],pch=15,col=colnum[j])}
    legend('topright',paste('alpha=',round(alpha,2),sep=""),lty=1,,lwd=3,col=colnum)
  }

  return(list(cv.lambda=min.lambda, cv.alpha=min.alpha,cv.mse=mseest[minloc],mse=mseest, lambda=lambdaest, alpha=alphaest,z=z))

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
