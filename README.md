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
source("data_gen.R", encoding = "UTF-8")

### generate data
# exactly sparse when approx_sparse=FALSE, approximately sparse when approx_sparse=TRUE
# setting 1
data_list = gen_data(n=1000,p=1500,setting="1",approx_sparse=FALSE)

# setting 3 and 4
# data_list = gen_data(n=1000,p=1500,setting="3")
# data_list = gen_data(n=1000,p=1500,setting="4",df=10) # df:degree of freedom for t distribution

X = data_list$X
y = data_list$y

### true value
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

```


### Setting 2
```R
library(DLL)
source("data_gen.R", encoding = "UTF-8")

### generate data
data_list = gen_data(n=1000,p=1500,setting="2")
X = data_list$X
y = data_list$y

### true value
f.deriv <- function(x) {
  mean.x = -0.25
  sd.x = 1
  p.x = function(x) pnorm(x,mean = mean.x, sd = sd.x)
  d.x = function(x) dnorm(x,mean = mean.x, sd = sd.x)
  return(-1.5*pi*cos(pi*p.x(x))*d.x(x))
}
f.deriv(d0)

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# D.ind=1: inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

```


### Comparison with ReSmoothing Estimator
The source code to implement ReSmoothing estimator is contained in the file spaddinf.R, and this is a direct copy of the original implementation by the author of this method https://github.com/gregorkb/spaddinf/tree/master/R.
```R
library(DLL)
library(grplasso)
source("data_gen.R", encoding = "UTF-8")
source("spaddinf.R", encoding = "UTF-8")

### generate data
# exactly sparse when approx_sparse=FALSE, approximately sparse when approx_sparse=TRUE
data_list = gen_data(n=500,p=750,setting="1", approx_sparse=FALSE)
X = data_list$X
y = data_list$y

### true value
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI
```


### Nonlinear Treatment Model
```R
library(DLL)
source("data_gen.R", encoding = "UTF-8")

### generate data
data_list = gen_data(n=1000,p=1500,setting="nonlinear_treatment")
X = data_list$X
y = data_list$y

### true value
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

### DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0, treatment.SAM = TRUE)

### point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI
```

