library(locpol)
library(nprobust)
library(SAM)

### plug-in local linear estimator
plug_in = function(X,y,D.ind,d0) {
  sam.model = cv.SAM(X,y,quant.trans=FALSE,lam.seq = NULL)
  f.hat = predict.SAM(sam.model,X)
  R.hat = y - apply(f.hat[,-(D.ind+1)], 1, sum) # skip intercept
  bw = suppressWarnings(thumbBw(X[,D.ind],R.hat,deg=1,kernel=SqK))
  # bw = lpbwselect(R.hat,X[,D.ind],eval=d0,deriv=1,kernel="uni",bwselect="mse-rot")$bws[,"h"]
  plug = lprobust(R.hat,X[,D.ind],eval=d0,deriv=1,p=1,h=bw,kernel="uni")
  est = plug$Estimate[,"tau.us"]
  est.se = plug$Estimate[,"se.us"]
  
  return(list(est = est, est.se = est.se))
}

### oracle local linear estimator
orac = function(X,y,e,g,D.ind,d0) {
  bw = suppressWarnings(thumbBw(X[,D.ind],g(X[,D.ind])+e,deg=1,kernel=SqK))
  # bw = lpbwselect(g(X[,D.ind])+e,X[,D.ind],eval=d0,deriv=1,kernel="uni",bwselect="mse-dpi")$bws[,"h"]
  orac = lprobust(g(X[,D.ind])+e,X[,D.ind],eval=d0,deriv=1,p=1,h=bw,kernel="uni")
  est = orac$Estimate[,"tau.us"]
  est.se = orac$Estimate[,"se.us"]
  
  return(list(est = est, est.se = est.se))
}
