#' Conduct sparse mediation for parge p ( p > n) with L1 or L2 penalization using fast computation of inverse matrix
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
#' V=500
#' set.seed(1234)
#'     covmat=matrix(0,V+2,V+2);
#'     covmat[1,2]=0.5;covmat[1, (1:3)+2]=rep(r.a,3);covmat[2, (1:3)+2]=rep(r.b,3);
#'     covmat=covmat+t(covmat);diag(covmat)<-1
#'     sqrtmat = sqrtmat.comp(covmat)
#'     tmpmat = matrix(rnorm(N*(V+2)),N,V+2) %*% sqrtmat
#'     X=tmpmat[,1]
#'     Y=tmpmat[,2]
#'     M=tmpmat[,-c(1:2)]
#'     sparse.mediation.largep(X,M,Y,tol=10^(-10),max.iter=100,lambda = log(1+(1:25)/50))
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation L1penalization
#' @import parallel
#' @import MASS
#' @import glmnet
#' @export
sparse.mediation.largep = function(X,M,Y,tol=10^(-10),max.iter=10,lambda = log(1+(1:50)/125),
                                glmnet.penalty.factor=c(0,rep(1,2*V)),alpha=1,threshold=0.00001,verbose=FALSE){


  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
  N = nrow(M)
  V = ncol(M)
  #Y.mean=mean(Y)
  #X.mean=mean(X)
  #M.mean=apply(M,2,mean)
  Y.sd=as.vector(sqrt(var(Y)))
  X.sd=as.vector(sqrt(var(X)))
  M.sd=sqrt(apply(M,2,var))
  Y = scale(Y,center=TRUE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=TRUE,scale=TRUE)

  ## Penalty Factor
  if (ncol(X)>1){stop("X has more than 1 colum. Stop.")}
  ## Initialization###
  ## OLS Estimation ###
  U = cbind(X,M)

  #invtMM = ginv(t(M)%*%M)
  tXX = t(X)%*%X
  tUY = t(U)%*%Y
  tMX = t(M)%*%X

  #tUU = #rbind(cbind(tXX, t(tMX)),cbind(tMX, t(M)%*%M))
  #tUU.sqmat=sqrtmat.comp(tUU)
  tUU = ginv.largep(U,sqrtmat=TRUE,sqrtinvmat=TRUE)


  ## Interative Update

  betaest =  matrix(0,1+2*V,length(lambda))

  for (j in 1:length(lambda)){
    if(verbose==TRUE){print(paste("Lambda",lambda[j]))}
    gamma_new = rep(0,V+1)#tUU$inv %*% tUY
    alpha_new = rep(0,V)#t(ginv(tXX)%*%t(X)%*%M)

    iter=0
    err=1000
    while( err>tol & iter<max.iter){

      alpha_old=alpha_new
      gamma_old = gamma_new
      beta_old = c(gamma_old,alpha_old)
#      print(beta_old)
      sigma1 = mean( (Y - U %*% gamma_old)^2)
      tmp = M - matrix(X,N,1) %*% matrix(alpha_old,1,V)
      #Sigma2 = t(tmp)%*%tmp/N

#      print('spot1')
      Sigma2=ginv.largep(tmp/sqrt(N),sqrtmat=TRUE,sqrtinvmat=TRUE)#ginv(Sigma2)

      #A = matrix(0,1+2*V,1+2*V)
      #A[1:(1+V),1:(1+V)]=1/sigma1 * tUU
      #A[(1+V)+ 1:V,(1+V)+ 1:V]=  as.numeric(tXX) * Sigma2.inv
      Asqmat = bdiag(1/sqrt(sigma1) * tUU$sqrtmat,  sqrt(as.numeric(tXX)) * Sigma2$sqrtinv)
      Asqmat.inv=bdiag(sqrt(sigma1) * tUU$sqrtinv,  1/sqrt(as.numeric(tXX)) * Sigma2$sqrtmat)

      #print(Asqmat)
      #print(Asqmat.inv)

      #Asqmat = A;Asqmat[1:(1+V),1:(1+V)]=1/sqrt(sigma1) * tUU.sqmat
      #Asqmat[(1+V)+ 1:V,(1+V)+ 1:V]=  sqrt(as.numeric(tXX)) * Sigma2.sqrt.inv
      #C = ginv(Asqmat) %*% rbind(tUY/sigma1, Sigma2.inv%*%tMX)
      C = Asqmat.inv %*% rbind(tUY/sigma1, Sigma2$inv%*%tMX)

      if(is.null(glmnet.penalty.factor)==TRUE){
        fit = glmnet(as.matrix(Asqmat), as.matrix(C),lambda=lambda[j],alpha=alpha)
      }else{
        fit = glmnet(as.matrix(Asqmat), as.matrix(C),lambda=lambda[j],penalty.factor=glmnet.penalty.factor,alpha=alpha)
      }

      beta_new=coef(fit)[-1]

      if (threshold>0){beta_new[abs(beta_new)<threshold]<-0}

      ## use thresholds as well: since all variables are standardized, coefficients less than 0.001 does not have any meaning.

      #beta_new[(1:V) +1]*beta_new[(1:V) +V+1]
      gamma_new = beta_new[1:(V+1)]
      alpha_new = beta_new[(1:V)+ V+1]
      err = sum((beta_old[-1]-beta_new[-1])^2)
      iter=iter+1
      if (verbose==TRUE){print(c(iter, err))}
    }
    betaest[,j]=beta_new
  }
  cest =betaest[1,]
  medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]
  nump=apply(betaest,2,function(x){sum(abs(x)>0)})


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
