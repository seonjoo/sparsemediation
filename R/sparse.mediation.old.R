#' Conduct sparse medication with elastic net (Old version)
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
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:50)/125)) tuning parameter for L1 penalization
#' @param alpha (defult=1) tuning parameter for L2 penalization
#' @param figure (defult=NULL) print figures for mean predictive errors by tuning parameters alpha and lambda
#' @param glmnet.penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @return c directeffect
#' @return hatb Path b (M->Y given X) estimates
#' @return hata Path a (X->M) estimates
#' @return medest Mediation estimates (a*b)
#' @return alpha
#' @return lambda
#' @return nump Number of selected mediation paths
#' @examples
#' N=100
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' sparse.mediation.old(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:25)/50))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @export

sparse.mediation.old = function(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:50)/125),
                            figure=NULL,display=FALSE,
                            glmnet.penalty.factor=c(0,rep(1,2*V)),alpha=1){
  library(parallel)
  library(MASS)
  library(glmnet)
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

  #if (display==TRUE){
  #  cols=rainbow(V)
    #if(is.null(figure)==FALSE){png(paste(figure,'.png',sep=""),width=1600,height=1600,res=200)}
  #  par(mfrow=c(2,2))
  #  x1 =  betaest[(1:V)+1,];x2=betaest[(1:V)+V+1,]
  #  plot(lambda,x1[1,],type='l',ylim=range(x1),lty=1,col=cols[1], main='Beta',ylab='beta',xlab='Lambda Loc')
  #  for (j in 2:V){lines(lambda,x1[j,],lty=1,col=cols[j])}
  #  plot(lambda,x2[1,],type='l',ylim=range(x2),lty=1,col=cols[1], main='Alpha',ylab='alpha',xlab='Lambda Loc')
  #  for (j in 2:V){lines(lambda,x2[j,],lty=2,col=cols[j])}
  #  plot(lambda,medest[1,],type='l',ylim=range(medest),lty=1,col=cols[1],main='Mediation',ylab='alpha*beta',xlab='Lambda Loc')
  #  for (j in 2:V){lines((lambda),medest[j,],lty=1,col=cols[j])}

  #  plot(lambda,nump,pch=15,main='Number of Non-zero Parameters');lines(lambda,nump)
    #if(is.null(figure)==FALSE){dev.off()}
  #  }

  return(list(
    c = cest,
    hatb=betaest[(1:V)+1,]*Y.sd/M.sd,
    hata=betaest[(1:V)+V+1,]*M.sd/X.sd,
    medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]*Y.sd/X.sd,
    alpha=alpha,
    lambda = lambda,
    nump=nump
  ))
}
