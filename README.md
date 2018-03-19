# sparsemediation
Conducting sparse medition

## From github:
```{r}
library(devtools)
install_github("seonjoo/sparsemediation")
library(sparsemediation)
```

## Usage
The sparse.mediation function conducts sparse mediation for the specified tuning parameters for elastic net. 
Cross-validation is possible using cv.sparse.mediation.

```{r}
N=100
V=50
set.seed(1234)
a = rbinom(V,1,0.1)*5;b<-a
X = scale(rnorm(N))
M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
Y =  10*X + M %*% b + rnorm(N)

cvfit<-cv.sparse.mediation(X, M, Y, tol = 10^(-10), K = 4, max.iter = 100,lambda = log(1 + (1:10)/25), tau=c(0.5,1,2),  multicore = 4, seednum = 1e+06)

fit<-sparse.mediation(X,M,Y,tol=10^(-10),max.iter=100,lambda = cvfit$cv.lambda,tau=cvfit$cv.tau)
