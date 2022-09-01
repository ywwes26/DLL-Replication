# DLL-Replication
This repository provides codes and guidance to replicate the simulation results of Decorrelated Local Linear estimator proposed in \<arxiv:1907.12732\>. For replication of the results, it is also helpful to check the package github page https://github.com/zijguo/HighDim-Additive-Inference and [Reference Manual](https://cran.r-project.org/web/packages/DLL/DLL.pdf).

## Installation of the Package DLL
The pakcage DLL can be installed from [CRAN](https://cran.r-project.org/):
```R
install.packages("DLL")
```

## Simulation Settings and Replication

### Setting 1, 3 and 4
```R
library(DLL)
source("gen_data.R", encoding = "UTF-8")
source("helpers.R", encoding = "UTF-8") # sparse additive model
source("local_linear.R", encoding = "UTF-8")

### generate data
# exactly sparse when approx_sparse=FALSE, approximately sparse when approx_sparse=TRUE
# setting 1
data_list = gen_data(n=1500,p=1500,setting="1",approx_sparse=FALSE)

# setting 3 and 4
# data_list = gen_data(n=1000,p=1500,setting="3")
# data_list = gen_data(n=1000,p=1500,setting="4",df=10) # df: degree of freedom for t distribution

X = data_list$X
y = data_list$y
e = data_list$e

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# index of D in the matrix(D,X)
D.ind = 1
# inference on the D.ind component of X
DLL.out = DLL(X=X, y=y, D.ind=D.ind, d0=d0)

### true value, changes as index of D changes, check gen_data.R
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

### plug-in local lienar estimator
plug.out = plug_in(X,y,D.ind,d0)
# point estimates and se
plug.out$est
plug.out$est.se

### oracle local linear estimator
orac.out = orac(X,y,e,g=function(x) 1.5*sin(x),D.ind,d0)
# point estimates and se
orac.out$est
orac.out$est.se
```


### Setting 2
```R
library(DLL)
source("gen_data.R", encoding = "UTF-8")

### generate data
data_list = gen_data(n=1000,p=1500,setting="2")
X = data_list$X
y = data_list$y

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# index of D in the matrix(D,X)
D.ind = 1
# inference on the D.ind component of X
DLL.out = DLL(X=X, y=y, D.ind=D.ind, d0=d0)

### true value, changes as index of D changes, check gen_data.R
f.deriv <- function(x) {
  mean.x = -0.25
  sd.x = 1
  p.x = function(x) pnorm(x,mean = mean.x, sd = sd.x)
  d.x = function(x) dnorm(x,mean = mean.x, sd = sd.x)
  return(-1.5*pi*cos(pi*p.x(x))*d.x(x))
}
f.deriv(d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

```


### Comparison with ReSmoothing Estimator
The source code to implement ReSmoothing estimator is contained in the file spaddinf.R, and this is a direct copy of the original implementation by the author of this method https://github.com/gregorkb/spaddinf/tree/master/R. This file is copied here since some internal functions not outputed in the package need to be used.
```R
library(DLL)
library(grplasso)
library(nprobust)
source("gen_data.R", encoding = "UTF-8")
source("spaddinf.R", encoding = "UTF-8")

### generate data
# exactly sparse when approx_sparse=FALSE, approximately sparse when approx_sparse=TRUE
data_list = gen_data(n=500,p=750,setting="1", approx_sparse=FALSE)
X = data_list$X
y = data_list$y

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# index of D in the matrix(D,X)
D.ind = 1
# inference on the D.ind component of X
DLL.out = DLL(X=X, y=y, D.ind=D.ind, d0=d0)

### true value, changes as index of D changes, check gen_data.R
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

### ReSmoothing estimator
# d.pre=20 is the default choice of the method
spaddinf.presmt.cv.out = spadd.presmth.Bspl.cv(X,y,d.pre=20,n.lambda=25,n.eta=25,n.folds=5)
# inference on the 1:D.ind components
spaddinf.presmt.out = spadd.presmth.Bspl(X,y,d.pre=20,lambda=spaddinf.presmt.cv.out$cv.lambda,eta=spaddinf.presmt.cv.out$cv.eta,n.foi=D.ind)

# the same bandwidth selection method as DLL
# the D.ind column of f.hat.design is the presmoothing estimator of the D.ind component
bw = suppressWarnings(thumbBw(X[,D.ind],spaddinf.presmt.out$f.hat.design[,D.ind],deg=1,kernel=SqK))
# use local linear on the presmoothing estimator
rs = lprobust(spaddinf.presmt.out$f.hat.design[,D.ind],X[,D.ind],eval=d0,deriv=1,p=1,h=bw,kernel="uni")
# point estimates and se
rs$Estimate[,"tau.us"]
rs$Estimate[,"se.us"]
```


### Nonlinear Treatment Model
```R
library(DLL)
source("gen_data.R", encoding = "UTF-8")

### generate data
data_list = gen_data(n=1000,p=1500,setting="nonlinear_treatment")
X = data_list$X
y = data_list$y

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0, treatment.SAM = TRUE)

### true value, changes as index of D changes, check gen_data.R
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI
```

