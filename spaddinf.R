#' Helps build the basis functions for Legendre polynomials
pcws.poly <- function(x,left.endpoints,K){
  
  vec <- numeric(length(left.endpoints)*(K+1))
  int <- sum(x >= left.endpoints) # in which interval does it fall?
  vec[((int-1)*(K+1)+1):((int)*(K+1))] <- 1*(x - left.endpoints[int])^c(0:K)
  return(vec)
  
}

#' Fit the desparsified lasso presmoothing estimator with cubic B-splines
#'
#' @param Y the response vector (centered)
#' @param X the design matrix
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X}) for which to compute the desparsified lasso presmoothing estimator.
#' @return a list with the fitted functions etc.
#' @examples
#' data <- data_gen(n = 200, q = 10, r = .5)
#'
#' spadd.presmth.Bspl.out <- spadd.presmth.Bspl(X = data$X,
#'                                              Y = data$Y,
#'                                              d.pre = 20,
#'                                              lambda = 1,
#'                                              eta = 3,
#'                                              n.foi = 6)
#'
#' plot_presmth_Bspl(spadd.presmth.Bspl.out,
#'                   true.functions = list( f = data$f,
#'                                          X = data$X))
#' @export
spadd.presmth.Bspl <- function(X,Y,d.pre,lambda,eta,n.foi)
{
  
  q <- ncol(X)
  n <- nrow(X)
  
  # make cubic B-splines basis function design
  HH <- HH.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length=q)
  knots.list <- vector("list",length=q)
  emp.cent <- vector("list",length=q)
  
  for( j in 1:q )
  {
    
    int.knots <- quantile(X[,j],seq(0,1,length=d.pre-2+1)) # add one, so that one can be removed after centering to restore full-rank.
    boundary.knots <- range(int.knots)
    all.knots <- sort(c(rep(boundary.knots,3),int.knots))
    knots.list[[j]] <- all.knots
    
    Bj <- spline.des(all.knots,X[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
    emp.cent[[j]] <- apply(Bj,2,mean)
    Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre,byrow=TRUE)
    
    # construct matrix in which l2 norm of function is a quadratic form
    M <- t(Bj.cent) %*% Bj.cent / n
    
    Q <- chol(M)
    Q.inv <- solve(Q)
    QQ.inv[[j]] <- Q.inv
    
    # construct basis function matrices
    HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
    HH <- cbind(HH,Bj.cent)
    groups <- c(groups,rep(j,d.pre))
    
  }
  
  # get the group lasso estimators
  grplasso.out <- grplasso(	y = Y,
                            x = HH.tilde,
                            index = groups,
                            lambda = lambda,
                            model = LinReg(),
                            center = FALSE,
                            standardize = FALSE,
                            control = grpl.control(trace=0))
  
  beta.tilde <- grplasso.out$coef
  lasso.fitted.values <- grplasso.out$fitted
  selected <- grplasso.out$norms.pen != 0
  
  Y.lasso <- f.hat.design <- matrix(NA,nrow(X),n.foi)
  f.hat <- vector("list",n.foi)
  AAt <- vector("list",n.foi)
  sigma.hat <- numeric(n.foi)
  
  for(j in 1:n.foi)
  {
    
    # Now get Z1, the matrix replacing the projection of H_1 onto the
    # orthogonal complement of the other columns of H.  Use the group
    # lasso to produce a projection having orthogonality
    # controlled by lambda.
    ind <- which(groups == j)
    
    Zj <- matrix(0,n,d.pre)
    
    for(l in 1:d.pre)
    {
      Z.model <- grplasso(y = HH[, ind[l]],
                          x = HH.tilde[,-ind],
                          index = groups[-ind],
                          lambda = eta,
                          model = LinReg(),
                          center = FALSE,
                          standardize = FALSE,
                          control = grpl.control( trace = 0 )
      )
      
      Zj[,l] <- HH[,ind[l]] - Z.model$fitted
      
    }
    
    # construct the presmoothing estimator
    Xj <- HH[,ind]
    Yj.lasso <- Y - HH.tilde[,-ind] %*% beta.tilde[-ind]
    
    Aj <- Xj %*% solve( t(Zj) %*% Xj ) %*% t(Zj)
    AAt[[j]] <- Aj %*% t(Aj)
    
    betaj.hat <- solve( t(Zj) %*% Xj ) %*% t(Zj) %*% Yj.lasso
    fj.hat <- Xj %*% betaj.hat	# presmoothing estimator
    
    f.hat.design[,j] <- fj.hat
    Y.lasso[,j] <- Yj.lasso
    
    # export actual function estimate
    f.hat[[j]] <- eval(parse(text=paste("function(x)
                                        {
                                        x <- round(x,10)
                                        x.mat <- spline.des(",paste("c(",paste(round(knots.list[[j]],6),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                        x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent[[j]],collapse=","),"),length(x),",d.pre,sep=""),",byrow=TRUE)
                                        f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(betaj.hat,collapse=","),")",sep=""),")
                                        return(f.hat)
                                        }"
    )))
    
    # do variance estimation
    D <- cbind(diag(n-1),rep(0,n-1)) - cbind(rep(0,n-1),diag(n-1))
    
    Xj.sort <-  Xj[order(X[,j]),]
    Zj.sort <-  Zj[order(X[,j]),]
    Aj.sort <- Xj.sort %*% solve( t(Zj.sort) %*% Xj.sort) %*% t(Zj.sort)
    
    nu.hat <- sum(diag( t(D %*% Aj.sort ) %*% D %*% Aj.sort ))
    sigma.hat[j] <- sqrt(sum((D %*% fj.hat[order(X[,j])])^2) / nu.hat)
    
  }
  
  output <- list(	f.hat.design = f.hat.design,
                  f.hat = f.hat,
                  Y.lasso = Y.lasso,
                  knots.list = knots.list,
                  AAt = AAt,
                  sigma.hat = sigma.hat)
  
  return(output)
  
}

#' Fit the desparsified lasso presmoothing estimator with Legendre polynomials
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X}) for which to compute the desparsified lasso presmoothing estimator.
#' @param K the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.
#' @return a list with the fitted functions etc.
#'
#' @examples
#' data <- data_gen(n = 100,q = 10,r = .5)
#'
#' spadd.presmth.Legr.out <- spadd.presmth.Legr(X = data$X,
#'                                              Y = data$Y,
#'                                              d.pre = 10,
#'                                              lambda = .1,
#'                                              eta = 2,
#'                                              n.foi = 6,
#'                                              K = 1)
#'
#' plot_presmth_Legr(x = spadd.presmth.Legr.out,
#'                   true.functions = list( f = data$f, X = data$X))
#' @export
# spadd.presmth.Legr <- function(X,Y,d.pre,lambda,eta,n.foi,K=1)
# {
#   
#   q <- ncol(X)
#   n <- nrow(X)
#   
#   # make cubic B-splines basis function design
#   HH <- HH.tilde <- matrix(NA,n,0)
#   groups <- numeric()
#   QQ.inv <- vector("list",length=q)
#   left.endpoints.list <- vector("list",length=q)
#   emp.cent <- vector("list",length=q)
#   
#   for( j in 1:q )
#   {
#     
#     left.endpoints <- quantile(X[,j],seq(0,1,length=d.pre+2))[-(d.pre+2)] # add one interval so that one can be removed after centering to restore full-rank.
#     left.endpoints.list[[j]] <- left.endpoints
#     
#     # now make the Legendre polynomial basis
#     Bj <- t(sapply(X[,j],FUN=pcws.poly,left.endpoints = left.endpoints,K=K))[,-c(1:(K+1))] # remove columns corresponding to first interval
#     emp.cent[[j]] <- apply(Bj,2,mean)
#     Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre*(K+1),byrow=TRUE)
#     
#     # construct matrix in which l2 norm of function is a quadratic form
#     M <- t(Bj.cent) %*% Bj.cent / n
#     
#     Q <- chol(M)
#     Q.inv <- solve(Q)
#     QQ.inv[[j]] <- Q.inv
#     
#     # construct basis function matrices
#     HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
#     HH <- cbind(HH,Bj.cent)
#     groups <- c(groups,rep(j,d.pre*(K+1)))
#     
#   }
#   
#   # get the group lasso estimators
#   grplasso.out <- grplasso(	y = Y,
#                             x = HH.tilde,
#                             index = groups,
#                             lambda = lambda,
#                             model = LinReg(),
#                             center = FALSE,
#                             standardize = FALSE,
#                             control = grpl.control(trace=0))
#   
#   beta.tilde <- grplasso.out$coef
#   lasso.fitted.values <- grplasso.out$fitted
#   # selected <- grplasso.out$norms.pen != 0
#   
#   Y.lasso <- f.hat.design <- matrix(NA,nrow(X),n.foi)
#   f.hat <- vector("list",n.foi)
#   maxima <- numeric(n.foi)
#   AAt <- vector("list",n.foi)
#   sigma.hat <- numeric(n.foi)
#   
#   for(j in 1:n.foi)
#   {
#     
#     # Now get Z1, the matrix replacing the projection of H_1 onto the
#     # orthogonal complement of the other columns of H.  Use the group
#     # lasso to produce a projection having orthogonality
#     # controlled by lambda.
#     ind <- which(groups == j)
#     
#     Zj <- matrix(0,n,d.pre*(K+1))
#     
#     for(l in 1:(d.pre*(K+1)))
#     {
#       Z.model <- grplasso(y = HH[, ind[l]],
#                           x = HH.tilde[,-ind],
#                           index = groups[-ind],
#                           lambda = eta,
#                           model = LinReg(),
#                           center = FALSE,
#                           standardize = FALSE,
#                           control = grpl.control( trace = 0 )
#       )
#       
#       Zj[,l] <- HH[,ind[l]] - Z.model$fitted
#       
#     }
#     
#     
#     # construct presmoothing estimator
#     Xj <- HH[,ind]
#     Yj.lasso <- Y - HH.tilde[,-ind] %*% beta.tilde[-ind]
#     
#     Aj <- Xj %*% solve( t(Zj) %*% Xj ) %*% t(Zj)
#     AAt[[j]] <- Aj %*% t(Aj)
#     
#     betaj.hat <- solve( t(Zj) %*% Xj ) %*% t(Zj) %*% Yj.lasso
#     fj.hat <- Xj %*% betaj.hat	# presmoothing estimator
#     
#     f.hat.design[,j] <- fj.hat
#     Y.lasso[,j] <- Yj.lasso
#     
#     # export actual function estimate
#     f.hat[[j]] <- eval(parse(text=paste("function(x)
#                                         {
#                                         x <- round(x,10)
#                                         x.mat <- t(sapply(x,FUN=pcws.poly,left.endpoints = c(",paste(round(left.endpoints.list[[j]],10),collapse=","),"), K=",K,"))[,-c(1:(",K,"+1))]
#                                         x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent[[j]],collapse=","),"),length(x),",d.pre*(K+1),sep=""),",byrow=TRUE)
#                                         f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(betaj.hat,collapse=","),")",sep=""),")
#                                         return(f.hat)
#                                         }"
#     )))
#     
#     maxima[j] <- max(X[,j])
#     
#     # do variance estimation
#     D <- cbind(diag(n-1),rep(0,n-1)) - cbind(rep(0,n-1),diag(n-1))
#     
#     Xj.sort <-  Xj[order(X[,j]),]
#     Zj.sort <-  Zj[order(X[,j]),]
#     Aj.sort <- Xj.sort %*% solve( t(Zj.sort) %*% Xj.sort) %*% t(Zj.sort)
#     
#     nu.hat <- sum(diag( t(D %*% Aj.sort ) %*% D %*% Aj.sort ))
#     sigma.hat[j] <- sqrt(sum((D %*% fj.hat[order(X[,j])])^2) / nu.hat)
#     
#     
#   }
#   
#   output <- list(	f.hat.design = f.hat.design,
#                   f.hat = f.hat,
#                   Y.lasso = Y.lasso,
#                   left.endpoints.list = left.endpoints.list,
#                   maxima = maxima,
#                   AAt = AAt,
#                   sigma.hat = sigma.hat)
#   
#   return(output)
#   
# }

#' Choose tuning parameters of desparsified lasso presmoothing estimator with cubic B-splines
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param n.lambda the number of candidate lambda values
#' @param n.eta the number of candidate eta values
#' @param n.folds the number of crossvalidation folds
#' @return a list with the chosen values of the tuning parameters
#'
#' @examples
#' data <- data_gen(n = 200,q = 50,r = .9)
#'
#' spadd.presmth.Bspl.cv.out <- spadd.presmth.Bspl.cv(X = data$X,
#'                                                    Y = data$Y,
#'                                                    d.pre = 20,
#'                                                    n.lambda = 25,
#'                                                    n.eta = 25,
#'                                                    n.folds = 5,
#'                                                    plot = TRUE)
#'
#' spadd.presmth.Bspl.out <- spadd.presmth.Bspl(X = data$X,
#'                                              Y = data$Y,
#'                                              d.pre = 20,
#'                                              lambda = spadd.presmth.Bspl.cv.out$cv.lambda,
#'                                              eta = spadd.presmth.Bspl.cv.out$cv.eta,
#'                                              n.foi = 6)
#'
#' plot_presmth_Bspl(spadd.presmth.Bspl.out)
#' @export
spadd.presmth.Bspl.cv <- function(X,Y,d.pre,n.lambda,n.eta,n.folds,plot = FALSE)
{
  
  q <- ncol(X)
  n <- nrow(X)
  
  # make cubic B-splines basis function design
  HH <- HH.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length=q)
  knots.list <- vector("list",length=q)
  emp.cent <- vector("list",length=q)
  
  for( j in 1:q )
  {
    
    int.knots <- quantile(X[,j],seq(0,1,length=d.pre-2+1)) # add one, so that one can be removed after centering to restore full-rank.
    boundary.knots <- range(int.knots)
    all.knots <- sort(c(rep(boundary.knots,3),int.knots))
    knots.list[[j]] <- all.knots
    
    Bj <- spline.des(all.knots,X[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
    emp.cent[[j]] <- apply(Bj,2,mean)
    Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre,byrow=TRUE)
    
    # construct matrix in which l2 norm of function is a quadratic form
    M <- t(Bj.cent) %*% Bj.cent / n
    Q <- chol(M)
    Q.inv <- solve(Q)
    QQ.inv[[j]] <- Q.inv
    
    # construct basis function matrices
    HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
    HH <- cbind(HH,Bj.cent)
    groups <- c(groups,rep(j,d.pre))
    
  }
  
  lambda.max <- lambdamax(y = Y,
                          x = HH.tilde,
                          index = groups,
                          model = LinReg(),
                          center = FALSE,
                          standardize = FALSE,
                          control = grpl.control(trace=0))
  
  lambda.min <- .001 * lambda.max
  lambda.seq <- c(exp(log(lambda.min) + ((n.lambda+1):1)/(n.lambda+1) * ((log(lambda.max) - log(lambda.min)))))[-1]
  
  # create list of sets of indices indicating which observations are in each fold
  folds <- vector("list", n.folds)
  fold.size <- floor(n / n.folds)
  
  for(fold in 1:n.folds){
    
    folds[[fold]] <- ((fold-1)*fold.size + 1):(fold*fold.size)
    
  }
  
  if( floor(n / n.folds) != n/n.folds )
  {
    folds[[n.folds]] <- c(folds[[n.folds]],(fold*fold.size+1):n)
  }
  
  
  # do crossvalidation for lambda
  cvMSEP <- matrix(0,n.folds,n.lambda)
  
  for(fold in 1:n.folds)
  {
    
    fold.ind <- folds[[fold]]
    
    grplasso.out.fold <- grplasso(y = Y[-fold.ind],
                                  x = HH.tilde[-fold.ind,],
                                  index = groups,
                                  lambda = lambda.seq*(n.folds-1)/n.folds,
                                  model = LinReg(),
                                  center = FALSE,
                                  standardize = FALSE,
                                  control = grpl.control(trace=0))
    
    Y.fold.mat <- matrix(Y[fold.ind],length(fold.ind),n.lambda)
    Y.hat.fold.mat <- HH.tilde[fold.ind,] %*% grplasso.out.fold$coef
    
    cvMSEP[fold,] <- apply((Y.fold.mat - Y.hat.fold.mat)^2,2,mean)
    
    print(paste("lambda fold: ", fold ,sep=""))
    
  }
  
  cv.MSEPs.lambda <- apply(cvMSEP,2,mean)
  which.lambda <- which.min(cv.MSEPs.lambda)
  cv.lambda <- lambda.seq[which.lambda]
  
  # crossvalidation for choosing eta
  # Do only for j = 1
  
  j <- 1
  ind <- which(groups == j)
  eta.max.vals <- numeric()
  for(l in 1:d.pre){
    
    eta.max.vals[l] <- lambdamax(y = HH[,ind[l]],
                                 x = HH.tilde[,-ind],
                                 index = groups[-ind],
                                 model = LinReg(),
                                 center = FALSE,
                                 standardize = FALSE,
                                 control = grpl.control(trace=0))
  }
  
  eta.max <- max(eta.max.vals)
  
  eta.min <- .001 * eta.max
  eta.seq <- c(exp(log(eta.min) + ((n.eta+1):1)/(n.eta+1) * ((log(eta.max) - log(eta.min)))))[-1]
  
  
  PRED.eta <- array(0,dim=c(n,d.pre,n.eta))
  cv.MSEP.eta <- array(0,dim=c(n.folds,n.eta,d.pre))
  cv.MSEP.eta.avg <- matrix(0,n.folds,n.eta)
  
  for(fold in 1:n.folds)
  {
    
    fold.ind <- folds[[fold]]
    
    for(l in 1:d.pre)
    {
      
      grplasso.out.fold <- grplasso( y = HH[-fold.ind,ind[l]],
                                     x = HH.tilde[-fold.ind,-ind],
                                     index = groups[-ind],
                                     lambda = eta.seq*(n.folds-1)/n.folds,
                                     model = LinReg(),
                                     center = FALSE,
                                     standardize = FALSE,
                                     control = grpl.control(trace=0))
      
      Y.fold.mat <- matrix(HH[fold.ind,ind[l]],length(fold.ind),n.eta)
      Y.hat.fold.mat <- HH.tilde[fold.ind,-ind] %*% grplasso.out.fold$coef
      
      cv.MSEP.eta[fold,,l] <- apply((Y.fold.mat - Y.hat.fold.mat)^2,2,mean)
      
    }
    
    cv.MSEP.eta.avg[fold,] <- apply(cv.MSEP.eta[fold,,],1,mean)
    
    print(paste("eta fold: ", fold ,sep=""))
    
  }
  
  cv.MSEPs.eta <- apply(cv.MSEP.eta.avg,2,mean)
  which.eta <- which.min(cv.MSEPs.eta)
  cv.eta <- eta.seq[which.eta]
  
  if(plot == TRUE)
  {
    
    # par(mfrow=c(1,2))
    
    plot(cv.MSEPs.lambda~log(lambda.seq),ylim=range(cv.MSEPs.lambda))
    abline(v = log(cv.lambda))
    
    plot(cv.MSEPs.eta~log(eta.seq),ylim=range(cv.MSEPs.eta))
    abline(v = log(cv.eta))
    
  }
  
  output <- list( cv.lambda = cv.lambda,
                  cv.eta = cv.eta)
  
  return(output)
  
}


#' Choose tuning parameters for fitting the desparsified lasso presmoothing estimator with Legendre polynomials
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param n.lambda the number of candidate lambda values
#' @param n.eta the number of candidate eta values
#' @param n.folds the number of crossvalidation folds
#' @param K the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.
#' @param plota logical indicating whether crossvalidation output should be plotted
#' @return a list with the chosen values of the tuning parameters
#'
#' @examples
#' data <- data_gen(n = 200,q = 50,r = .9)
#'
#' spadd.presmth.Legr.cv.out <- spadd.presmth.Legr.cv(X = data$X,
#'                                                    Y = data$Y,
#'                                                    d.pre = 10,
#'                                                    n.lambda = 25,
#'                                                    n.eta = 25,
#'                                                    n.folds = 5,
#'                                                    plot = TRUE)
#'
#' spadd.presmth.Legr.out <- spadd.presmth.Legr(X = data$X,
#'                                              Y = data$Y,
#'                                              d.pre = 10,
#'                                              lambda = spadd.presmth.Legr.cv.out$cv.lambda,
#'                                              eta = spadd.presmth.Legr.cv.out$cv.eta,
#'                                              n.foi = 6)
#'
#' plot_presmth_Legr(x = spadd.presmth.Legr.out,
#'                   true.functions = list( f = data$f,
#'                                          X = data$X))
#' @export
# spadd.presmth.Legr.cv <- function(X,Y,d.pre,K = 1,n.lambda,n.eta,n.folds,plot = FALSE)
# {
#   
#   q <- ncol(X)
#   n <- nrow(X)
#   
#   # make cubic B-splines basis function design
#   HH <- HH.tilde <- matrix(NA,n,0)
#   groups <- numeric()
#   QQ.inv <- vector("list",length=q)
#   left.endpoints.list <- vector("list",length=q)
#   emp.cent <- vector("list",length=q)
#   
#   for( j in 1:q )
#   {
#     
#     left.endpoints <- quantile(X[,j],seq(0,1,length=d.pre+2))[-(d.pre+2)] # add one interval so that one can be removed after centering to restore full-rank.
#     left.endpoints.list[[j]] <- left.endpoints
#     
#     # now make the Legendre polynomial basis
#     Bj <- t(sapply(X[,j],FUN=pcws.poly,left.endpoints = left.endpoints,K=K))[,-c(1:(K+1))] # remove columns corresponding to first interval
#     emp.cent[[j]] <- apply(Bj,2,mean)
#     Bj.cent <- Bj - matrix(emp.cent[[j]],n,d.pre*(K+1),byrow=TRUE)
#     
#     # construct matrix in which l2 norm of function is a quadratic form
#     M <- t(Bj.cent) %*% Bj.cent / n
#     
#     Q <- chol(M)
#     Q.inv <- solve(Q)
#     QQ.inv[[j]] <- Q.inv
#     
#     # construct basis function matrices
#     HH.tilde <- cbind(HH.tilde,Bj.cent %*% Q.inv)
#     HH <- cbind(HH,Bj.cent)
#     groups <- c(groups,rep(j,d.pre*(K+1)))
#     
#   }
#   
#   lambda.max <- lambdamax(y = Y,
#                           x = HH.tilde,
#                           index = groups,
#                           model = LinReg(),
#                           center = FALSE,
#                           standardize = FALSE,
#                           control = grpl.control(trace=0))
#   
#   lambda.min <- .001 * lambda.max
#   lambda.seq <- c(exp(log(lambda.min) + ((n.lambda+1):1)/(n.lambda+1) * ((log(lambda.max) - log(lambda.min)))))[-1]
#   
#   # create list of sets of indices indicating which observations are in each fold
#   folds <- vector("list", n.folds)
#   fold.size <- floor(n / n.folds)
#   
#   for(fold in 1:n.folds){
#     
#     folds[[fold]] <- ((fold-1)*fold.size + 1):(fold*fold.size)
#     
#   }
#   
#   if( floor(n / n.folds) != n/n.folds )
#   {
#     folds[[n.folds]] <- c(folds[[n.folds]],(fold*fold.size+1):n)
#   }
#   
#   
#   # do crossvalidation for lambda
#   cvMSEP <- matrix(0,n.folds,n.lambda)
#   
#   for(fold in 1:n.folds)
#   {
#     
#     fold.ind <- folds[[fold]]
#     
#     grplasso.out.fold <- grplasso(y = Y[-fold.ind],
#                                   x = HH.tilde[-fold.ind,],
#                                   index = groups,
#                                   lambda = lambda.seq*(n.folds-1)/n.folds,
#                                   model = LinReg(),
#                                   center = FALSE,
#                                   standardize = FALSE,
#                                   control = grpl.control(trace=0))
#     
#     Y.fold.mat <- matrix(Y[fold.ind],length(fold.ind),n.lambda)
#     Y.hat.fold.mat <- HH.tilde[fold.ind,] %*% grplasso.out.fold$coef
#     
#     cvMSEP[fold,] <- apply((Y.fold.mat - Y.hat.fold.mat)^2,2,mean)
#     
#     print(paste("lambda fold: ", fold ,sep=""))
#     
#   }
#   
#   cv.MSEPs.lambda <- apply(cvMSEP,2,mean)
#   which.lambda <- which.min(cv.MSEPs.lambda)
#   cv.lambda <- lambda.seq[which.lambda]
#   
#   # crossvalidation for choosing eta
#   # Do only for j = 1
#   
#   j <- 1
#   ind <- which(groups == j)
#   eta.max.vals <- numeric()
#   for(l in 1:(d.pre*(K+1))){
#     
#     eta.max.vals[l] <- lambdamax(y = HH[,ind[l]],
#                                  x = HH.tilde[,-ind],
#                                  index = groups[-ind],
#                                  model = LinReg(),
#                                  center = FALSE,
#                                  standardize = FALSE,
#                                  control = grpl.control(trace=0))
#   }
#   
#   eta.max <- max(eta.max.vals)
#   
#   eta.min <- .001 * eta.max
#   eta.seq <- c(exp(log(eta.min) + ((n.eta+1):1)/(n.eta+1) * ((log(eta.max) - log(eta.min)))))[-1]
#   
#   PRED.eta <- array(0,dim=c(n,d.pre*(K+1),n.eta))
#   cv.MSEP.eta <- array(0,dim=c(n.folds,n.eta,d.pre*(K+1)))
#   cv.MSEP.eta.avg <- matrix(0,n.folds,n.eta)
#   
#   for(fold in 1:n.folds)
#   {
#     
#     fold.ind <- folds[[fold]]
#     
#     for(l in 1:(d.pre*(K+1)))
#     {
#       
#       grplasso.out.fold <- grplasso( y = HH[-fold.ind,ind[l]],
#                                      x = HH.tilde[-fold.ind,-ind],
#                                      index = groups[-ind],
#                                      lambda = eta.seq*(n.folds-1)/n.folds,
#                                      model = LinReg(),
#                                      center = FALSE,
#                                      standardize = FALSE,
#                                      control = grpl.control(trace=0))
#       
#       Y.fold.mat <- matrix(HH[fold.ind,ind[l]],length(fold.ind),n.eta)
#       Y.hat.fold.mat <- HH.tilde[fold.ind,-ind] %*% grplasso.out.fold$coef
#       
#       cv.MSEP.eta[fold,,l] <- apply((Y.fold.mat - Y.hat.fold.mat)^2,2,mean)
#       
#     }
#     
#     cv.MSEP.eta.avg[fold,] <- apply(cv.MSEP.eta[fold,,],1,mean)
#     
#     print(paste("eta fold: ", fold ,sep=""))
#     
#   }
#   
#   cv.MSEPs.eta <- apply(cv.MSEP.eta.avg,2,mean)
#   which.eta <- which.min(cv.MSEPs.eta)
#   cv.eta <- eta.seq[which.eta]
#   
#   if(plot == TRUE)
#   {
#     
#     # par(mfrow=c(1,2))
#     
#     plot(cv.MSEPs.lambda~log(lambda.seq),ylim=range(cv.MSEPs.lambda))
#     abline(v = log(cv.lambda))
#     
#     plot(cv.MSEPs.eta~log(eta.seq),ylim=range(cv.MSEPs.eta))
#     abline(v = log(cv.eta))
#     
#   }
#   
#   output <- list( cv.lambda = cv.lambda,
#                   cv.eta = cv.eta)
#   
#   return(output)
#   
# }


#' Fit simple nonparametric regression model with cubic B-splines
#'
#' @param Y a response vector
#' @param X vector of covariate observations
#' @param d the number functions in the cubic B-spline basis
#' @return a list containing the fitted function and a vector containing the values of the fitted function at the design points
#'
#' @examples
#' data <- data_gen(n = 200, q = 50, r = .9)
#'
#' spadd.presmth.Bspl.out <- spadd.presmth.Bspl(X = data$X,
#'                                              Y = data$Y,
#'                                              d.pre = 20,
#'                                              lambda = 1,
#'                                              eta = 3,
#'                                              n.foi = 6)
#'
#' resmth.Bspl.out <- resmth.Bspl(Y = data$Y.oracle[,1],
#'                                X = data$X[,1],
#'                                d = 6,
#'                                AAt = spadd.presmth.Bspl.out$AAt[[1]],
#'                                sigma.hat = spadd.presmth.Bspl.out$sigma.hat[1],
#'                                plot = TRUE)
#' @export
resmth.Bspl <- function(Y,X,d,AAt,sigma.hat,plot = FALSE,x = NULL,derivs=0,alpha = 0.05)
{
  
  if( length(x) == 0 ) x <- seq(min(X)+1e-2,max(X)-1e-2,length = 200)
  
  n <- length(Y)
  
  int.knots <- quantile(X,seq(0,1,length=d-2+1)) # add one, so that one can be removed after centering to restore full-rank.
  boundary.knots <- range(int.knots)
  all.knots <- sort(c(rep(boundary.knots,3),int.knots))
  
  B <- spline.des(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
  emp.cent <- apply(B,2,mean)
  B.cent <- B - matrix(emp.cent,n,d,byrow=TRUE)
  
  BB.cent.inv <- solve( t(B.cent) %*% B.cent)
  beta.hat <- as.numeric(BB.cent.inv %*% t(B.cent) %*% (Y - mean(Y)))
  f.hat.design <- as.numeric(B.cent %*% beta.hat)
  
  x.mat <- spline.des(all.knots, x = x, ord = 4, derivs = derivs, outer.ok = TRUE)$design[,-1]
  x.mat.cent <- x.mat - matrix(emp.cent,length(x),d,byrow=TRUE)
  f.hat.x <- x.mat.cent %*% beta.hat
  
  V <- BB.cent.inv %*% t(B.cent) %*% AAt %*% B.cent %*% BB.cent.inv
  f.hat.se.x <- sqrt(diag(x.mat.cent %*% V %*% t(x.mat.cent)))
  
  CIu.x <- f.hat.x + qnorm(1-alpha/2) * f.hat.se.x * sigma.hat
  CIl.x <- f.hat.x - qnorm(1-alpha/2) * f.hat.se.x * sigma.hat
  
  if(plot == TRUE)
  {
    
    plot(Y ~ X, col = rgb(0.545,0,0,1))
    lines(f.hat.x ~ x,lwd=1.5,col=rgb(0,0,.545))
    
    x.poly <- c(x,x[length(x):1])
    y.poly <- c(CIl.x,CIu.x[length(x):1])
    polygon(x = x.poly, y = y.poly, col = rgb(0,0,.545,.5),border=NA)
    
    
  }
  
  output <- list(f.hat.design = f.hat.design,
                 f.hat.x = f.hat.x,
                 CIu.x = CIu.x,
                 CIl.x = CIl.x,
                 d = d,
                 sigma.hat = sigma.hat,
                 alpha = alpha)
  
  return(output)
  
}


#' Choose number of basis functions in least squares B-splines via crossvalidation
#'
#' @param Y a response vector (centered)
#' @param X vector of covariate observations
#' @param d.seq a sequence of candidate numbers of basis functions
#' @param n.folds the number of crossvalidation folds
#' @param plot a logical indicating whether to plot crossvalidation output
#' @return the number of basis functions selected by crossvalidation
#' @examples
#' data <- data_gen(n = 200,
#'                  q = 50,
#'                  r = .9)
#'
#' cv.d <- Bspl.cv(Y = data$Y.oracle[,1],
#'                 X = data$X[,1],
#'                 d.seq = 3:12,
#'                 n.folds = 5,
#'                 plot = TRUE)
#'
#' oracle.Bspl.out <- oracle.Bspl(Y = data$Y.oracle[,1],
#'                                X = data$X[,1],
#'                                d = cv.d,
#'                                plot = TRUE)
#' @export
Bspl.cv  <- function(Y,X,d.seq, n.folds, plot = FALSE)
{
  
  n <- length(X)
  cv.msep <- numeric()
  
  # create list of sets of indices indicating which observations are in each fold
  folds <- vector("list", n.folds)
  fold.size <- floor(n / n.folds)
  
  # reorder indices so that we do not remove a bunch of contiguous observations
  ind.reordered <- as.numeric((matrix(order(X),fold.size,n.folds,byrow=TRUE)))
  
  for(fold in 1:n.folds){
    
    folds[[fold]] <- ind.reordered[((fold-1)*fold.size + 1):(fold*fold.size)]
    
  }
  
  if( floor(n / n.folds) != n/n.folds )
  {
    folds[[n.folds]] <- ind.reordered[c(folds[[n.folds]],(fold*fold.size+1):n)]
  }
  
  # go through different numbers of basis functions
  cvMSEP <- matrix(0,n.folds,length(d.seq))
  for(l in 1:length(d.seq))
  {
    
    int.knots <- quantile(X,seq(0,1,length=d.seq[l]-2+1)) # add one, so that one can be removed after centering to restore full-rank.
    boundary.knots <- range(int.knots)
    all.knots <- sort(c(rep(boundary.knots,3),int.knots))
    
    B <- spline.des(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
    emp.cent <- apply(B,2,mean)
    B.cent <- B - matrix(emp.cent,n,d.seq[l],byrow=TRUE)
    
    Y.hat <- numeric(n)
    
    # go through crossvalidation folds
    for( fold in 1:n.folds)
    {
      
      fold.ind <- folds[[fold]]
      
      Y.hat[fold.ind] <- B.cent[fold.ind,] %*% solve( t(B.cent[-fold.ind,]) %*% B.cent[-fold.ind,]) %*% t(B.cent[-fold.ind,]) %*% Y[-fold.ind]
      
      cvMSEP[fold,l] <- mean( (Y[fold.ind] - Y.hat[fold.ind])^2 )
      
    }
    
  }
  
  cv.MSEPs <- apply(cvMSEP,2,mean)
  which.d <- which.min(cv.MSEPs)
  cv.d <- d.seq[which.d]
  
  
  output <- list(cv.d)
  
  
  if(plot == TRUE)
  {
    
    plot(cv.MSEPs ~ d.seq)
    
  }
  
  return(cv.d)
  
}


#' Fit two step estimator with B-splines in the first step and B-splines in the second step
#'
#' @param Y the response vector (centered)
#' @param X the design matrix
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param d.re the number of intervals in which to divide the support of each covariate for the resmoother
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X}) for which to compute the desparsified lasso presmoothing estimator.
#' @param x a sequence of values at which the final estimators should be evaluated
#' @param K the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.
#' @return a list with the fitted functions and pointwise confidence intervals
#'
#' @examples
#' data <- data_gen(n = 200, q = 50, r = .9)
#'
#' preresmth.Bspl.Bspl.out <- preresmth.Bspl.Bspl(Y = data$Y,
#'                                                X = data$X,
#'                                                d.pre = 20,
#'                                                d.re = 10,
#'                                                lambda = 5,
#'                                                eta = 3,
#'                                                n.foi = 6,
#'                                                alpha = 0.05)
#' @export
preresmth.Bspl.Bspl <- function(X,Y,d.pre,d.re = NULL,lambda,eta,n.foi,x.eval,derivs=0,Y.oracle,plot=FALSE,alpha = 0.05)
{
  
  spadd.presmth.Bspl.out <- spadd.presmth.Bspl(X = X,
                                               Y = Y,
                                               d.pre = d.pre,
                                               lambda = lambda,
                                               eta = eta,
                                               n.foi = n.foi)
  
  f.hat.x <- CIl.x <- CIu.x <- f.hat.x.ori <- CIl.x.ori <- CIu.x.ori <- matrix(0,length(x.eval),n.foi)
  cv.d <- numeric(n.foi)
  
  # if(length(d.re) == 0){
  #   
  #   
  #   d[j] <- Bspl.cv(Y = spadd.presmth.Bspl.out$f.hat.design[,j],
  #                   X = X[,j],
  #                   d.seq = 5:floor(d.pre*7/8),
  #                   n.folds = 5,
  #                   plot = FALSE)
  #   
  # } else {
  #   
  #   d <- rep(d.re,n.foi)
  #   
  # }
  
  
  for(j in 1:n.foi)
  {
    
    cv.d[j] <- Bspl.cv(Y = Y.oracle[,j],
                       X = X[,j],
                       d.seq = 5:floor(d.pre*7/8),
                       n.folds = 5,
                       plot = FALSE)
    
    
    resmth.Bspl.out <- resmth.Bspl(Y = Y.oracle[,j],
                                   X = X[,j],
                                   d = cv.d[j],
                                   AAt = spadd.presmth.Bspl.out$AAt[[j]],
                                   sigma.hat = spadd.presmth.Bspl.out$sigma.hat[j],
                                   plot = plot,
                                   x = x.eval,
                                   derivs = 1,
                                   alpha = alpha)
    
    resmth.Bspl.ori.out <- resmth.Bspl(Y = Y.oracle[,j],
                                   X = X[,j],
                                   d = cv.d[j],
                                   AAt = spadd.presmth.Bspl.out$AAt[[j]],
                                   sigma.hat = spadd.presmth.Bspl.out$sigma.hat[j],
                                   plot = plot,
                                   x = x.eval,
                                   derivs = 0,
                                   alpha = alpha)
    
    f.hat.x[,j] <- resmth.Bspl.out$f.hat.x
    CIl.x[,j] <- resmth.Bspl.out$CIl.x
    CIu.x[,j] <- resmth.Bspl.out$CIu.x
    f.hat.x.ori[,j] <- resmth.Bspl.ori.out$f.hat.x
    CIl.x.ori[,j] <- resmth.Bspl.ori.out$CIl.x
    CIu.x.ori[,j] <- resmth.Bspl.ori.out$CIu.x
    
  }
  
  output <- list( f.hat.x = f.hat.x,
                  CIl.x = CIl.x,
                  CIu.x = CIu.x,
                  f.hat.x.ori = f.hat.x.ori,
                  CIl.x.ori = CIl.x.ori,
                  CIu.x.ori = CIu.x.ori,
                  x = x.eval,
                  sigma.hat = spadd.presmth.Bspl.out$sigma.hat,
                  cv.d = cv.d,
                  alpha = alpha
  )
  
  return(output)
  
}


#' Fit two step estimator with Legendre polynomials in first step and B-splines in the second step
#'
#' @param Y the response vector (centered)
#' @param X the design matrix
#' @param d.pre the number of intervals in which to divide the support of each covariate
#' @param lambda the tuning parameter for fitting the group lasso estimate for the bias correction
#' @param eta the tuning parameter for the group lasso projection of one set of basis functions onto those of the other covariates.
#' @param n.foi the number of functions (first columns of \code{X}) for which to compute the desparsified lasso presmoothing estimator.
#' @param x a sequence of values at which the final estimators should be evaluated
#' @param K the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.
#' @return a list with the fitted functions and pointwise confidence intervals
#'
#' @examples
#' data <- data_gen(n = 200, q = 50, r = .9)
#'
#' preresmth.Legr.Bspl.out <- preresmth.Legr.Bspl(Y = data$Y,
#'                                                X = data$X,
#'                                                d.pre = 20,
#'                                                d.re = 10,
#'                                                lambda = 5,
#'                                                eta = 3,
#'                                                n.foi = 6,
#'                                                plot = TRUE,
#'                                                alpha = 0.05)
#' @export
# preresmth.Legr.Bspl <- function(X,Y,d.pre,d.re = NULL,lambda,eta,n.foi,plot=FALSE,alpha = 0.05)
# {
#   
#   spadd.presmth.Legr.out <- spadd.presmth.Legr(X = X,
#                                                Y = Y,
#                                                d.pre = d.pre,
#                                                lambda = lambda,
#                                                eta = eta,
#                                                n.foi = n.foi)
#   
#   f.hat.x <- CIl.x <- CIu.x <- matrix(0,200,n.foi)
#   cv.d <- numeric(n.foi)
#   xx <- matrix(0,200,n.foi)
#   
#   if(length(d.re) == 0){
#     
#     
#     d[j] <- Bspl.cv(Y = spadd.presmth.Legr.out$f.hat.design[,j],
#                     X = X[,j],
#                     d.seq = 5:floor(d.pre*7/8),
#                     n.folds = 5,
#                     plot = FALSE)
#     
#   } else {
#     
#     d <- rep(d.re,n.foi)
#     
#   }
#   
#   
#   for(j in 1:n.foi)
#   {
#     
#     cv.d[j] <- Bspl.cv(Y = spadd.presmth.Legr.out$f.hat.design[,j],
#                        X = X[,j],
#                        d.seq = 5:floor(d.pre*7/8),
#                        n.folds = 5,
#                        plot = FALSE)
#     
#     xx[,j] <- seq(min(X[,j]),max(X[,j]),length = 200)
#     
#     resmth.Bspl.out <- resmth.Bspl(Y = spadd.presmth.Legr.out$f.hat.design[,j],
#                                    X = X[,j],
#                                    d = d[j],
#                                    AAt = spadd.presmth.Legr.out$AAt[[j]],
#                                    sigma.hat = spadd.presmth.Legr.out$sigma.hat[j],
#                                    plot = plot,
#                                    x = xx[,j],
#                                    alpha = alpha)
#     
#     f.hat.x[,j] <- resmth.Bspl.out$f.hat.x
#     CIl.x[,j] <- resmth.Bspl.out$CIl.x
#     CIu.x[,j] <- resmth.Bspl.out$CIu.x
#     
#   }
#   
#   output <- list( f.hat.x = f.hat.x,
#                   CIl.x = CIl.x,
#                   CIu.x = CIu.x,
#                   x = xx,
#                   sigma.hat = spadd.presmth.Legr.out$sigma.hat,
#                   d = d,
#                   alpha = alpha
#   )
#   
#   return(output)
#   
# }

#' Fit simple nonparametric regression model with cubic B-splines
#'
#' @param Y a response vector (centered)
#' @param X vector of covariate observations
#' @param d the number functions in the cubic B-spline basis
#' @param a logical indicating whether a plot of the fitted function should be generated
#' @param x a vector of values at which evaluations of the fitted function should be returned/plotted
#' @return a list containing the fitted function and a vector containing the values of the fitted function at the design points
#'
#' @examples
#' data <- data_gen(n = 200, q = 50, r = .9)
#'
#' oracle.Bspl.out <- oracle.Bspl(Y = data$Y.oracle[,1],
#'                                X = data$X[,1],
#'                                d = 10,
#'                                plot = TRUE,
#'                                x = NULL)
#' @export
oracle.Bspl <- function(Y,X,d,plot = FALSE,x = NULL,derivs=0,alpha = 0.05)
{
  
  if( length(x) == 0 ) x <- seq(min(X)+1e-2,max(X)-1e-2,length = 200)
  
  n <- length(Y)
  
  int.knots <- quantile(X,seq(0,1,length=d-2+1)) # add one, so that one can be removed after centering to restore full-rank.
  boundary.knots <- range(int.knots)
  all.knots <- sort(c(rep(boundary.knots,3),int.knots))
  
  B <- spline.des(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
  emp.cent <- apply(B,2,mean)
  B.cent <- B - matrix(emp.cent,n,d,byrow=TRUE)
  
  BB.cent.inv <- solve( t(B.cent) %*% B.cent)
  beta.hat <- as.numeric(BB.cent.inv %*% t(B.cent) %*% (Y - mean(Y)))
  f.hat.design <- as.numeric(B.cent %*% beta.hat)
  
  x.mat <- spline.des(all.knots,x = x, ord = 4, derivs = derivs, outer.ok = TRUE)$design[,-1]
  x.mat.cent <- x.mat - matrix(emp.cent,length(x),d,byrow=TRUE)
  f.hat.x <- x.mat.cent %*% beta.hat
  
  f.hat.se.x <- sqrt(diag(x.mat.cent %*% BB.cent.inv %*% t(x.mat.cent)))
  
  D <- cbind(diag(n-1),rep(0,n-1)) - cbind(rep(0,n-1),diag(n-1))
  sigma.hat <- sqrt(sum((D %*% Y[order(X)])^2) / (2*(n-1)))
  
  CIu.x <- f.hat.x + qnorm(1 - alpha/2) * f.hat.se.x * sigma.hat
  CIl.x <- f.hat.x - qnorm(1 - alpha/2) * f.hat.se.x * sigma.hat
  
  if(plot == TRUE)
  {
    
    plot(Y ~ X)
    lines(f.hat.x ~ x,lwd=1.5,col=rgb(0,0,.545))
    
    x.poly <- c(x,x[length(x):1])
    y.poly <- c(CIl.x,CIu.x[length(x):1])
    polygon(x = x.poly, y = y.poly, col = rgb(0,0,.545,.5),border=NA)
    
    
  }
  
  output <- list(f.hat.design = f.hat.design,
                 f.hat.x = f.hat.x,
                 CIu.x = CIu.x,
                 CIl.x = CIl.x,
                 d = d,
                 sigma.hat = sigma.hat,
                 alpha = alpha)
  
  return(output)
}

#' Generate data for simulation studies
#' @export
data_gen <- function(n,q,r)
{
  
  f <- vector("list",q)
  f[[1]] <- function(x){-sin(x*2)}
  f[[2]] <- function(x){x^2 - 25/12}
  f[[3]] <- function(x){x}
  f[[4]] <- function(x){exp(-x)-2/5*sinh(5/2)}
  f[[5]] <- function(x){x*0}
  f[[6]] <- function(x){x*0}
  
  # get design matrix
  R <- r^abs( outer(1:q,1:q,"-"))
  P <- 2*sin( R * pi / 6)
  X <- (pnorm( matrix(rnorm(n*q),ncol = q) %*% chol(P)) - .5) * 5
  
  signal.uncent <- cbind(f[[1]](X[,1]),
                         f[[2]](X[,2]),
                         f[[3]](X[,3]),
                         f[[4]](X[,4]),
                         f[[5]](X[,5]),
                         f[[6]](X[,6]))
  
  signal.means <- apply(signal.uncent,2,mean)
  signal <- signal.uncent - matrix(signal.means,n,6,byrow=TRUE)
  noise <- rnorm(n)
  Y <- apply(signal,1,sum) + noise - mean(noise)
  
  Y.oracle <- signal + matrix(noise - mean(noise),n,6,byrow=FALSE)
  
  output <- list(X = X,
                 Y = Y,
                 Y.oracle = Y.oracle,
                 signal.means = signal.means,
                 f = f,
                 n = n,
                 q = q,
                 r = r)
  
  return(output)
  
}


#' Plot the first six fitted functions
#' @export
plot_presmth_Bspl <- function(x, true.functions = NULL)
{
  
  f.hat.design <- x$f.hat.design
  f.hat <- x$f.hat
  knots.list <- x$knots.list
  
  par(mfrow=c(2,3),mar=c(2.1,2.1,1.1,1.1))
  
  for( j in 1:6)
  {
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    
    plot(NA,
         ylim = range(f.hat.design),
         xlim=c(xj.min,xj.max),
         ylab="",
         xlab="")
    
    plot(f.hat[[j]],min(xj.min),max(xj.max),col=rgb(.545,0,0,1),add=TRUE)
    
    if(length(true.functions)!=0)
    {
      
      x.seq <- seq(xj.min,xj.max,length=300)
      f.cent.seq <- true.functions$f[[j]](x.seq) - mean(true.functions$f[[j]](true.functions$X[,j]))
      lines(f.cent.seq ~ x.seq,lty=2)
    }
    
  }
  
}


#' plot presmoother when it is made of Legendre polynomials
#' @export
plot_presmth_Legr <- function(x,true.functions = NULL)
{
  
  f.hat <- x$f.hat
  f.hat.design <- x$f.hat.design
  left.endpoints.list <- x$left.endpoints.list
  maxima <- x$maxima
  
  par(mfrow=c(2,3),mar=c(2.1,2.1,1.1,1.1))
  
  for( j in 1:6)
  {
    
    
    xj.min <- min(left.endpoints.list[[j]]) + 1e-2
    xj.max <- max(left.endpoints.list[[j]]) - 1e-2
    
    
    plot(NA,
         xlim=c(xj.min,maxima[j]),
         ylim=range(f.hat.design),
         xlab="",
         ylab="")
    
    # points(f.hat.design[order(X[,j]),j]~sort(X[,j]))
    # plot each fitted Legendre polynomial function
    for(l in 1:(length(left.endpoints.list[[j]])-1))
    {
      
      x1 <- left.endpoints.list[[j]][l]
      x2 <- left.endpoints.list[[j]][l+1]
      plot(f.hat[[j]],x1,x2-1e-4,col=rgb(.545,0,0,1),lwd=1.5,add=TRUE)
    }
    
    plot(f.hat[[j]],x2,maxima[j]-1e-2,col=rgb(.545,0,0,1),lwd=1.5,add=TRUE)
    
    if(length(true.functions)!=0)
    {
      
      x.seq <- seq(xj.min,xj.max,length=300)
      f.cent.seq <- true.functions$f[[j]](x.seq) - mean(true.functions$f[[j]](true.functions$X[,j]))
      lines(f.cent.seq ~ x.seq,lty=2)
    }
    
  }
}