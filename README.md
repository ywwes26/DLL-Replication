# DLL-Replication
This repository provides codes and guidance to replicate the simulation results of Decorrelated Local Linear estimator proposed in \<arxiv:1907.12732\>.

## Installation of the Package DLL
The pakcage DLL can be installed from [CRAN](https://cran.r-project.org/):
```R
install.packages("DLL")
```

## Simulation Settings and Replication

### Setting 1
```R
library(DLL)
source("data_gen.R", encoding = "UTF-8")

# generate data
data_list = data_gen(n=1000,p=1500,setting="1",approx_sparse=FALSE)
data_list = data_gen(n=1000,p=1500,setting="1",approx_sparse=FALSE)
X = data_list$X
y = data_list$y

# DLL estimator
# evaluation points
d0 = c(0.1, 0.25)
# inference on the first component of X
DLL.out = DLL(X=X, y=y, D.ind=1, d0=d0)

# true value
f.deriv = function(d) 1.5*cos(d)
f.deriv(d0)

# point estimates, se and CI
DLL.out$est
DLL.out$est.se
DLL.out$CI

```
