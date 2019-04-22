### the support function file

###################### 10/18/2017 modify the support function: getXY since the information matrix


######################################################################################################
############################### the function to calculate the CV2 value ##############################
######################################################################################################
### this function is discarded.
cv.value.lsa <- function(sh.resp,CV=5,iter=5,low=0.1,upp=0.9,V = F,Lambda = 10^seq(-5,3,by=0.1), cvf = 1, Initial = NULL){
  cv <- cv.sd <- c(); p <- ncol(oz); n <- nrow(oz)
  for(i in 1:iter){
    mgs <- cv.lsa2(sh.resp,CV=CV,low=low, upp=upp, V = V, cv = cvf, Initial = Initial, loss = loss)
    cv <- rbind(cv, mgs$Loss.V)
  }
  Loss.sum <- apply(cv,2,sum)/nrow(cv)
  Loss.sd <- apply(cv,2,sd)/sqrt(nrow(cv))
  min.ind <- which(Loss.sum == min(Loss.sum))
  cv.min <- Lambda[min.ind]
  se.ind <- which(Loss.sum >= (Loss.sum[min.ind]+Loss.sd[min.ind]))[1]
  cv.1se <- Lambda[se.ind]
  
  list(loss = loss, cv.sd = cv.sd, cv.min = cv.min, cv.1se = cv.1se)
}


####################re-arrange the covariance matrix as the design matrix#######################
### see section 3.3, step 2 and step 3
getVY <- function(cov,n,coef,tm){######construct the design matrix
  p <- ncol(coef)
  if(is.null(p)) p<-1
  
  tmd <- tm[-1] - head(tm, -1) # since the time difference at last time point is 0
  tmd <- c(tmd, median(tmd)) # since the time difference at last time point is 0, we use median of time diff to represent it
  
  V.inv <- lapply(1:n, function(i) solve(cov[-1,-1,i]))
  #V.inv <- lapply(1:n, function(i) solve(cov[,,i]))
  V.inv <- Matrix::bdiag(V.inv) 
  ###transforme the V into the order for each path, see Appendix A
  eigen.V <- eigen(V.inv)
  sigma.root <- eigen.V$vectors %*% sqrt(diag(eigen.V$values)) %*% t(eigen.V$vectors)
  design <- sqrt(diag(rep(tmd,each = p))) %*% sigma.root
  #V.w <- diag(tmd) %*% eigen.V$vectors %*% sqrt(solve(diag(eigen.V$values))) %*% solve(eigen.V$vectors) #do not consider \delta(t)
  order <- c()
  for(o in 1:p){
    order <- c(order,seq(1,p*n,by=p)+o-1)
  }
  design <- design[order,order]
  #design <- design[-c(1:n), -c(1:n)] #since we doesn't consider intercept
  y.w <- design %*% as.vector(coef)
  list(Design = design, Y = y.w, V = V.inv)
}


#########################the function for penalty matrix(consider this part later) 09/03/2017 ########################
betad.fun <- function(tm, coef, lb, ub){
  n <- dim(coef)[1]; p <- dim(coef)[2]
  if(is.null(n)){
    n <- length(coef)
    p <- 1
    
    betad <- coef[-1] - head(coef, -1)
    tmd <- tm[-1] - head(tm, -1)
    
    sup <- max((betad/tmd)[lb:ub])#####exclude the unstable parts
    W1 <- diag(1/tmd) ###1/timedifference
    W2 <- diag(sum(sup)/sup/tmd) ###1/sup(first derivative)
    W3 <- diag(1/betad) ###1/estimated first derivative
    W4 <- diag(sum(sup)/(betad*sup))
    
    D1 <- W1 %*% getD1d(length(tm))
    D2 <- W2 %*% getD1d(length(tm))
    D3 <- W3 %*% getD1d(length(tm))
    D4 <- W4 %*% getD1d(length(tm))
  }
  else{
    betad <- coef[-1] - head(coef, -1) # a matrix
    tmd <- tm[-1] - head(tm, -1) # a vector
    
    sup <- apply(betad[lb:ub,], 2, max)
    W1 <- diag(1/tmd) ###1/timedifference
    W2 <- diag(sum(sup)/sup[1]/tmd,length(tmd)) ###1/sup(first derivative)
    W3 <- diag(1/betad[,1]) ###1/estimated first derivative
    W4 <- diag(sum(sup)/(betad[,1]*sup[1]*tmd))
    
    D1 <- W1 %*% getD1d(length(tm))
    D2 <- W2 %*% getD1d(length(tm))
    D3 <- W3 %*% getD1d(length(tm))
    D4 <- W4 %*% getD1d(length(tm))
    
    for(j in 2:p){
      W1 <- diag(1/tmd) ###1/timedifference
      W2 <- diag(sum(sup)/sup[j]/tmd,length(tmd)) ###1/sup(first derivative)
      W3 <- diag(1/betad[,j]) ###1/estimated first derivative
      W4 <- diag(sum(sup)/(betad[,j]*sup[j]*tmd))
      
      D1 <- bdiag(D1, W1 %*% getD1d(length(tm)))
      D2 <- bdiag(D2, W2 %*% getD1d(length(tm)))
      D3 <- bdiag(D3, W3 %*% getD1d(length(tm)))
      D4 <- bdiag(D4, W4 %*% getD1d(length(tm)))
    }
  }
  
  list(D1=D1, D2=D2, D3=D3, D4=D4, W1=W1, W2=W2, W3=W3, W4=W4, betad=betad, tmd=tmd)
}

#########fit the solution path for different interval##########
### see section 3.3, step 4
fit.genlasso <- function(mod, Design = F){
  n <- nrow(mod$oz);p <- ncol(mod$oz)
  #lower <- round(quantile(1:n,low));upper <- round(quantile(1:n,upp))
  risk <-  which(mod$unique.del==1)
  #lb <- sum(mod$delta[1:lower]); ub <- sum(mod$delta[1:upper])
  
  D = sd(mod$oz[,2]) * getD1d(length(risk))
  if(p > 2){#have more than one covariate
    for(i in 3:p){
      D <- bdiag(D,sd(mod$oz[,i]) * getD1d(length(risk)))
    }
  }
  D <- as.matrix(D)
  
  
  if(Design == T){
    Design.C <- getVY(mod$covmatrix[,,risk],length(risk),mod$coef[risk,-1], mod$time[risk])
    Y <- Design.C$Y
    X <- Design.C$Design
    #W.C <- betad.fun(mod$time[risk.c], mod$coef[risk.c,-1],lb=lb,ub=ub)
    
    out.gen <- genlasso(y=Y, X=X, D = D)###get the smoothed path for each lambda
    #out.gen2 <- genlasso(y=Y.C, X=X.C, D=as.matrix(W.C$D1))
    #out.gen3 <- genlasso(y=Y.C, X=X.C, D=as.matrix(W.C$D2))
    #out.gen4 <- genlasso(y=Y.C, X=X.C, D=as.matrix(W.C$D3))
    #out.gen5 <- genlasso(y=Y.C, X=X.C, D=as.matrix(W.C$D4))
    list(out = out.gen, risk=risk)
    #list(out.unwe = out.gen, out.td=out.gen2, out.sup=out.gen3,out.we=out.gen4, out.mix=out.gen5,
    #     range=c(low:upp),risk.c=risk.c)
  }
  else{
    Y <- as.vector(mod$coef[risk,-1])
    X <- diag(length(Y))
    W <- betad.fun(mod$time[risk], mod$coef[risk,-1],lb=lb,ub=ub)
    
    out.gen <- genlasso(y=Y, X=X, D=D)###get the smoothed path for each lambda
    out.gen2 <- genlasso(y=Y, X=X, D=as.matrix(W$D1))
    out.gen3 <- genlasso(y=Y, X=X, D=as.matrix(W$D2))
    out.gen4 <- genlasso(y=Y, X=X, D=as.matrix(W$D3))
    out.gen5 <- genlasso(y=Y, X=X, D=as.matrix(W$D4))
    list(out.unwe = out.gen, out.td=out.gen2, out.sup=out.gen3,out.we=out.gen4, out.mix=out.gen5,
         range=c(low:upp),risk=risk)
  }
}

### see section 3.3, step 5
restore.intercept <- function(sh.resp, lambda){
  p <- ncol(sh.resp$coef) -1
  cov <- sh.resp$covmatrix[,,sh.resp$unique.del==1]
  coef.unpe <- sh.resp$coef[sh.resp$unique.del==1,]
  whole.path <- fit.genlasso(sh.resp, Design =T)
  
  path.beta <- matrix(coef(whole.path$out, lambda = lambda)$beta, ncol = p)
  information <- lapply(1:sum(sh.resp$unique.del), function(i) solve(cov[,,i])) 
  
  intercept.pe <- lapply(1:sum(sh.resp$unique.del), function(i) coef.unpe[i,1] - 1/information[[i]][1,1] * t(information[[i]][1,-1]) %*% 
                           (path.beta[i,] - coef.unpe[i,-1]))
  intercept <- unlist(intercept.pe)
  
  list(lambda = lambda, fused.path = cbind(intercept, path.beta), time = sh.resp$time[sh.resp$unique.del==1])
}