library(MASS)
library(mvtnorm)

gen_data = function(n=1000, p=1500, setting="1", approx_sparse = FALSE, df = NULL) {
  
  # dimension of (D,X)
  p_dx = p + 1
  
  # functions
  g1 = function(x) 1.5*sin(x)
  g2 = function(x) 2*exp(-x/2)
  g3 = function(x) (x-1)^2 - 25/12
  g4 = function(x) x - 1/3
  g5 = function(x) 0.75*x
  g6 = function(x) 0.5*x
  g7 = function(x) 0.4*x
  g8 = function(x) 0.3*x
  g9 = function(x) 0.2*x
  g10 = function(x) 0.1*x
  g11 = function(x) 0.1*sin(2*pi*x)
  g12 = function(x) 0.2*cos(2*pi*x)
  g13 = function(x) 0.3*(sin(2*pi*x))^2
  g14 = function(x) 0.4*(cos(2*pi*x))^3
  g15 = function(x) 0.5*(sin(2*pi*x))^3
  
  # covariance structure
  Cov_Matrix = toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p_dx-4)))
  
  # data generation
  if (setting == "1") {
    mean.x = -0.25
    X = mvrnorm(n,rep(mean.x,p_dx),Sigma = Cov_Matrix)
    e = rnorm(n,sd=1)
    # generate response Y
    if (approx_sparse==FALSE) {
      y = g1(X[,1]) + g2(X[,2]) + g3(X[,3]) + g4(X[,4]) + g5(X[,5]) + g6(X[,6]) + e
    } 
    if (approx_sparse==TRUE) {
      y = g1(X[,1]) + g2(X[,2]) + g3(X[,3]) + g4(X[,4]) + 
        g5(X[,5]) + g6(X[,6]) + g7(X[,7]) + g8(X[,8]) + 
        g9(X[,9]) + g10(X[,10]) + g11(X[,11]) + g12(X[,12]) + 
        g13(X[,13]) + g14(X[,14]) +g15(X[,15]) + e
      for (j in 16:p_dx) {
        y = y + (j-1)^(-1)*X[, j]
      }
    }
  }
  
  if (setting == "2") {
    g1 = function(x) -1.5*sin(pi*x)
    g2 = function(x) 2*exp(-x)
    g5 = function(x) x^3 - 1/2
    g6 = function(x) x/(1+x)
    
    mean.x = -0.25
    sd.x = 1
    # CDF of X
    p.x = function(x) pnorm(x,mean = mean.x, sd = sd.x)
    X = mvrnorm(n,rep(mean.x,p_dx),Sigma = Cov_Matrix)
    Z = apply(X,2,p.x)
    e = rnorm(n,sd=1)
    # generate response Y
    y = g1(Z[,1]) + g2(Z[,2]) + g3(Z[,3]) + g4(Z[,4]) + g5(Z[,5]) + g6(Z[,6]) + e
  }
  
  if (setting == "3") {
    mean.x = -0.25
    sd.x = 1
    # CDF of X
    p.x = function(x) pnorm(x,mean = mean.x, sd = sd.x)
    X = mvrnorm(n,rep(mean.x,p_dx),Sigma = Cov_Matrix)
    X = (apply(X,2,p.x)-0.5)*5 # create correlated Uniform X on (-2.5,2.5)
    e = rnorm(n,sd=1)
    # generate response Y
    y = g1(X[,1]) + g2(X[,2]) + g3(X[,3]) + g4(X[,4]) + g5(X[,5]) + g6(X[,6]) + e
  }
  
  if (setting == "4") {
    X = rmvt(n,Cov_Matrix*(df-2)/df,df=df)
    e = rnorm(n,sd=1)
    # generate response Y
    y = g1(X[,1]) + g2(X[,2]) + g3(X[,3]) + g4(X[,4]) + g5(X[,5]) + g6(X[,6]) + e
  }
  
  if (setting == "nonlinear_treatment") {
    mean.x = -0.25
    Cov_Matrix = toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p-4)))
    X = mvrnorm(n,rep(mean.x,p),Sigma = Cov_Matrix)
    D = -0.5*exp(-X[,1]/2) + 0.5*sin(X[,2]) + 0.25*X[,3]^2 - 0.5*X[,4] - 0.25*X[,5]^2 + 0.5*cos(X[,6]) - 0.25*exp(X[,7]/2) + 0.25*X[,8] + rnorm(n,sd=0.5)
    X = cbind(D,X)
    e = rnorm(n,sd=1)
    # generate response Y
    y = g1(X[,1]) + g2(X[,2]) + g3(X[,3]) + g4(X[,4]) + g5(X[,5]) + g6(X[,6]) + e
  }
  
  return(list(X=X,y=y,e=e))
  
}
