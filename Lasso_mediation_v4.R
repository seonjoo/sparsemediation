#######################################################################
#### Sparse Structural Equation for Mediation Analysis
####
#### Codes were originally writtey by Seonjoo Lee on 05 13 2015 ####
#### Additional contributers: Jingyan Wang (MPH in Biostatistics), Tested Codes, evaluate simulation runs
####
#### Update: 2015 12 30
####  1) Finalize codes and evaluation based on glmnet
####
######################################################################

require(glmnet)
require(MASS)
require(parallel)
sparse.mediation = function(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:50)/125),
		figure=NULL,display=FALSE,
		glmnet.penalty.factor=c(0,rep(1,2*V)),alpha=1){
  
  # fit a GLM with lasso or elasticnet regularization
  # Description
  #
  # Fit a mediation model via penalized maximum likelihood and structural equation model. 
  # The regularization path is computed for the lasso or elasticnet penalty at a grid of 
  # values for the regularization parameter lambda. Currently, mediation analysis is developed based on gaussian assumption.
  # Multiple Mediaton Model:
  #   (1) M = Xa + e1
  #   (2) Y = Xc' + Mb + e2
  # And in the optimization, we do not regularize c', due to the assumption of partial mediation. 
  #Usage 
  #
  #lasso.mediation = function(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:50)/125),
  #                           figure=NULL,display=FALSE,
  #                           glmnet.penalty.factor=c(0,rep(1,2*V)),alpha=1)
  #
  #Arguments
  # 
  # X : N-by-1 matrix
  
  
  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
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
	
	## Penalty Factor
	if (ncol(X)>1){stop("X has more than 1 colum. Stop.")}
	## Initialization###
	## OLS Estimation ###
	U = cbind(X,M)
	tUU = t(U)%*%U
	invtUU = ginv(tUU)
	invtMM = ginv(t(M)%*%M)
	tXX = t(X)%*%X
	tUY = t(U)%*%Y
	tMX = t(M)%*%X

	## Interative Update

	betaest =  matrix(0,1+2*V,length(lambda))
	for (j in 1:length(lambda)){
		#print(paste("Lambda",lambda[j]))
		gamma_old = gamma_new = invtUU %*% tUY
		sigma1 = mean((Y - U %*% gamma_old)^2)
		alpha_old = alpha_new = t(solve(t(X)%*%X)%*%t(X)%*%M)
		beta_old = beta_new = c(gamma_old,alpha_old)
		tmp = M - matrix(X,N,1) %*% matrix(alpha_old,1,V)
		Sigma2 = t(tmp)%*%tmp/N

		iter=0
		err=1000
		while( err>tol & iter<max.iter){
			beta_old = beta_new
			alpha_old=alpha_new
			gamma_old = gamma_new
			A = matrix(0,1+2*V,1+2*V)
			A[1:(1+V),1:(1+V)]=1/sigma1 * tUU
			A[(1+V)+ 1:V,(1+V)+ 1:V]=as.numeric(tXX) * ginv(Sigma2)

			tmp=svd(A)
			sqmatA = tmp$u %*% diag(sqrt(tmp$d)) %*% t(tmp$v)
			C = solve(sqmatA) %*% rbind(tUY/sigma1, ginv(Sigma2)%*%tMX)


			if(is.null(glmnet.penalty.factor)==TRUE){
				fit = glmnet(sqmatA, C,lambda=lambda[j],alpha=alpha)
			}else{ fit = glmnet(sqmatA, C,lambda=lambda[j],penalty.factor=glmnet.penalty.factor,alpha=alpha)
			}

			beta_new = as.vector(predict(fit,type="coef"))[-1]
			## use thresholds as well: since all variables are standardized, coefficients less than 0.001 does not have any meaning. 
			beta_new[abs(beta_new)<0.001]<-0 

			#beta_new[(1:V) +1]*beta_new[(1:V) +V+1]
			gamma_new = beta_new[1:(V+1)] 
			alpha_new = beta_new[(1:V)+ V+1]
			sigma1 = mean((Y - U %*% gamma_new)^2)
			tmp = M - matrix(X,N,1) %*% matrix(alpha_new,1,V)
			Sigma2 = t(tmp)%*%tmp/N


			err = sum((beta_old-beta_new)^2)
			iter=iter+1
		}
		betaest[,j]=beta_new
	}
  	cest =betaest[1,] 
	medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]
	nump=apply(betaest,2,function(x){sum(abs(x)>0)})
	
	if (display==TRUE){
	cols=rainbow(V)
	if(is.null(figure)==FALSE){png(paste(figure,'.png',sep=""),width=1600,height=1600,res=200)}
	par(mfrow=c(2,2))
	x1 =  betaest[(1:V)+1,];x2=betaest[(1:V)+V+1,]
	plot(lambda,x1[1,],type='l',ylim=range(x1),lty=1,col=cols[1], main='Beta',ylab='beta',xlab='Lambda Loc')
	for (j in 2:V){lines(lambda,x1[j,],lty=1,col=cols[j])}
	plot(lambda,x2[1,],type='l',ylim=range(x2),lty=1,col=cols[1], main='Alpha',ylab='alpha',xlab='Lambda Loc')
	for (j in 2:V){lines(lambda,x2[j,],lty=2,col=cols[j])}
	plot(lambda,medest[1,],type='l',ylim=range(medest),lty=1,col=cols[1],main='Mediation',ylab='alpha*beta',xlab='Lambda Loc')
	for (j in 2:V){lines((lambda),medest[j,],lty=1,col=cols[j])}
	
	plot(lambda,nump,pch=15,main='Number of Non-zero Parameters');lines(lambda,nump)
	if(is.null(figure)==FALSE){dev.off()}}

return(list(
            c = cest,
            hatb=betaest[(1:V)+1,], 
            hata=betaest[(1:V)+V+1,],
            medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,],
            alpha=alpha,
            lambda = lambda,
            nump=nump
))
}

## We abort idea of differentially weighting paths a and b in the mediation due to too arbitrary choice of parameters, and also standardization helps. 
sparse.mediation.twoway = function(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:30)/100),alpha=(0:4)/4,figure=NULL){
	re=new.env()
	re=as.list(re)
	for (j in 1:length(alpha)){
		re[[j]]=sparse.mediation(X,M,Y,tol=tol,max.iter=max.iter,lambda = lambda,alpha=alpha[j],figure=NULL,display=FALSE)
		re[[j]]$alpha=alpha[j]

	}
	return(re)
}

################################################
# K-fold Cross-Validation
################################################

cv.sparse.mediation= function(X,M,Y,tol=10^(-10),K=5,max.iter=100,lambda = log(1+(1:15)/50),alpha=(0:4)/4,figure=NULL,multicore=1,seednum=1000000,display=TRUE){
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


