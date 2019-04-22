#created at 09/02/2017
# use least square approximation itself as the CV-criteria function.

#modified at 09/03/2017 change the form of loss function: consider time diffeence at loss function and use D as penalty matrix
#modified at 09/04/2017 to accomodate tied failure time, test it tomorrow; combine the cv2 into cv1 function
#modified at 9/21/2017 to accomodate tied failure events; the modification at 09/04 is not good.
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/supportFunc.R"))


cv.lsa2 <- function(sh.resp,CV=5,low=0.05, upp=0.95, V = F, cv = 1, Initial = NULL, loss = 2){ # initial is the intiial value matrix(usually from cox model)
  #V=F means we use the design matrix from whole sample; V=T means we use the design matrix from training data
  #split the data
  oz <- sh.resp$oz
  oy <- sh.resp$oy
  odelta <- sh.resp$delta
  otpt <- sh.resp$time
  uni.delta <- sh.resp$unique.del
  
  n <- nrow(oz)
  p <- ncol(oz)
  CV <- CV
  Loss.V <- Loss.SD <- c()
  n.lower <- as.integer(quantile(1:n,low))
  n.upper <- as.integer(quantile(1:n,upp))
  s.q <- sample(n.lower:n.upper,n.upper-n.lower+1,replace=FALSE)
  n.tru <- length(s.q)
  count <- 0 #count the effective number of effective cv
  for (J in 1:CV) {
    cat("J = ",J,"\n")
    te.ind <- sort(s.q[(1+(J-1)*floor(n.tru/CV)):(J*floor(n.tru/CV))])#test data
    tr.ind <- sort(c(s.q[-((1+(J-1)*floor(n.tru/CV)):(J*floor(n.tru/CV)))],(1:(n.lower-1)),((n.upper+1):n)))#training
    #tr.ind <- sort(c(s.q[-((1+(J-1)*floor(n.tru/CV)):(J*floor(n.tru/CV)))],(1:(n.lower-1))))
    oz.tr <- oz[-te.ind,]
    oy.tr <- oy[-te.ind]
    odel.tr <- odelta[-te.ind]
    otpt.tr <- unique(oy.tr)
    tr.ind.o <- match(otpt.tr, otpt)
    # odel.tr.o <- uni.delta[tr.ind.o]
    n.tr <- nrow(oz.tr)
    n.tr.o <- length(otpt.tr)
    Initial.tr <- Initial[tr.ind.o,]
    
    oz.te <- oz[te.ind,]
    oy.te <- oy[te.ind]
    odel.te <- odelta[te.ind]
    otpt.te <- unique(oy.te)
    te.ind.o <- match(otpt.te, otpt)
    # odel.te.o <- uni.delta[te.ind.o]
    n.te <- nrow(oz.te)
    n.te.o <- length(otpt.te)
    Initial.te <- Initial[te.ind.o,]
    ###construct the index after deleting the tied time points, then works on the unique index
    
    n.o <- n.te.o + n.tr.o
    
    lower <- round(quantile(1:n.tr,low));upper <- round(quantile(1:n.tr,upp))
    lb <- sum(odel.tr[1:lower]); ub <- sum(odel.tr[1:upper])#lower bound and upper bound
    sh.tr <- cox_temporal_op(oz.tr,oy.tr,odel.tr,init.fun=bini.fun, B=0, initial = Initial.tr, loss = loss)
    sh.te <- cox_temporal_op(oz.te,oy.te,odel.te,init.fun=bini.fun, B=0, initial = Initial.te, loss = loss)
    
    odel.tr.o <- sh.tr$unique.del
    odel.te.o <- sh.te$unique.del
    
    if(sh.tr$conv==0){ #means that this fold of training data is solvable
      count <- count + 1
      if(V == T){
        sh.tr <- cox_temporal_sh(oz.tr,oy.tr,odel.tr,init.fun=bini.fun, B=200)
        cov.tr <- sh.tr$covmatrix
        coef.tr <- cbind(order = tr.ind.o, sh.tr$coef)
        coef.tr1 <- rbind(coef.tr,cbind(order = te.ind.o, matrix(0, n.te.o, p)))
        coef.tr1 <- coef.tr1[order(coef.tr1[,1]),]########the whole solution path including te.ind recording as 0 firstly
        ##########get the parameter estimate and covariance matrix at each time points for test data set
        coef.te <- matrix(numeric(0), n.te.o, p)
        cov.te <- array(numeric(0),dim = c(p,p,n.te.o))
        belong <- ((te.ind.o-1) %in% tr.ind.o)#check the belong to status, means if the test subject is the next one of training subject
        for(i in 1:n.te.o){
          if(belong[i]){
            coef.te[i,] <- coef.tr1[te.ind.o[i]-1,-1] #doesn't include ID 
            cov.te[,,i] <- cov.tr[,,which(tr.ind.o==te.ind.o[i]-1)]
          }
          else if(te.ind[i] != 1) { # for subjects not the first one as well as belong[i] = F
            coef.te[i,] <- coef.te[i-1,]
            cov.te[,,i] <- cov.te[,,i-1]
          }
          else {
            coef.te[i,] <- c(log(1/n),0,0)
            cov.te[,,i] <- cov.tr[,,1]
          }
        }
      }
      else { #use the design matrix from whole sample, no longer to resampling procedure
        cov.tr <- sh.resp$covmatrix[,,tr.ind.o]
        cov.te <- sh.resp$covmatrix[,,te.ind.o]
        coef.tr <- sh.tr$coef
        #coef.tr1 <- rbind(coef.tr,cbind(order = te.ind.o, matrix(0, n.te.o, p)))
        #coef.tr1 <- coef.tr1[order(coef.tr1[,1]),]########the whole solution path including te.ind
        ##########get the parameter estimate and covariance matrix at each time points for test data set
        coef.te <- matrix(numeric(0), n.te.o, p)
        belong1 <- ((te.ind.o-1) %in% tr.ind.o)#check if the test time is the one next the one in training set
        belong2 <- te.ind.o %in% tr.ind.o # check if training and test set share the same time point
        for(i in 1:n.te.o){
          if(belong2[i]){
            coef.te[i,] <- coef.tr[which(tr.ind.o == te.ind.o[i]),]
          }
          else if(belong1[i]){
            coef.te[i,] <- coef.tr[which(tr.ind.o == (te.ind.o[i]-1)),]
          }
          else if(te.ind.o[i] != 1) { # example: 5,6 in test set; 4, 7, 8, in training set
            coef.te[i,] <- coef.te[i-1,]
          }
          else {
            coef.te[i,] <- c(log(1/n),rep(0, p-1))
          }
        }
      }
      
      ###transformation of covariance matrix into genlasso form. See appendix A for details.
      ###here we only consider the uncensored time points
      Designmatrix <- getVY(cov.tr[,,sh.tr$unique.del==1],sum(sh.tr$unique.del),
                            sh.tr$coef[sh.tr$unique.del==1,-1], sh.tr$time[sh.tr$unique.del == 1])
      Y <- Designmatrix$Y # T^{1/2} V^{-1/2} \hat{\beta}
      X <- Designmatrix$Design #T^{1/2} V^{-1/2}
      #W <- betad.fun(sh.tr$time[sh.tr$delta==1], sh.tr$coef[sh.tr$delta==1,-1], lb=lb, ub=ub) #the weight matrix
      ########################################################################################
      
      D = sd(oz.tr[,2]) * getD1d(sum(sh.tr$unique.del)) #the sd of first covariate
      
      if(p > 2){#have more than one covariate, we need block diagonal l1 fused matrix
        for(i in 3:p){ # when p = 3, we have 2 covariates since the first covariate is the intercept
          D <- bdiag(D,sd(oz.tr[,i]) * getD1d(sum(sh.tr$unique.del)))
        }
      }
      D <- as.matrix(D)
      
      out.gen <- genlasso(y=Y, X=X, D=D) ### only the 1d fused matrix
      
      # if(weight==1) out.gen <- genlasso(y=Y, X=X, D=D) ### only the 1d fused matrix
      # else if(weight == 2) out.gen <- genlasso(y=Y, X=X, D=as.matrix(W$D1)) #
      # else if(weight == 3) out.gen <- genlasso(y=Y, X=X, D=as.matrix(W$D2))
      # else if(weight == 4) out.gen <- genlasso(y=Y, X=X, D=as.matrix(W$D3))
      # else if(weight == 5) out.gen <- genlasso(y=Y, X=X, D=as.matrix(W$D4))
      
      
      ##################################################################################################################
      ###################################### cross-validation part #####################################################
      Lambda <- 10^seq(-5,4,by=0.1)
      #Lambda <- seq(0, 200, by = 1)
      ##################the coef array of gen lasso for each lambda###########################
      coef.tr.lasso <- array(numeric(0),dim = c(n.tr.o,p,length(Lambda)))
      fail.ind.tr <- which(sh.tr$unique.del == 1) 
      censor.ind.tr <- which(sh.tr$unique.del == 0) 
      #order.v <- order(c(tr.ind.o[odel.tr.o==1], tr.ind.o[odel.tr.o==0],te.ind.o))
      for(i in 1:length(Lambda)){
        coef.tr.lasso[fail.ind.tr,,i] <- matrix(c(sh.tr$coef[sh.tr$unique.del==1,1],coef(out.gen, lambda = Lambda[i])$beta), 
                                          sum(odel.tr.o), p)
        #coef.v.lasso[,,i] <- coef.v.lasso[order.v,,i]
      }
     for(i in censor.ind.tr){
       if(i == 1) coef.tr.lasso[1,,] <- c(log(1/n.tr), rep(0, p-1))
       else coef.tr.lasso[i,,] <- coef.tr.lasso[i-1,,]
     }
      
      ##################complete the solution path for censored subjects in training data#################
      # for(i in sort(c(tr.ind.o[odel.tr.o==0], te.ind.o))){
      #   if(i == 1) coef.v.lasso[i,,] <- rep(c(log(1/n),0,0), length(Lambda))
      #   else coef.v.lasso[i,,] <- coef.v.lasso[i-1,,]
      # }
      ##################the coef array for corresponding testing data set for each lambda#####
      coef.te.lasso <- array(numeric(0),dim = c(n.te.o,p,length(Lambda)))
      #coef.te.lasso <- coef.v.lasso[sort(te.ind.o),,]
      for(i in 1:n.te.o){
        if(belong2[i]){
          coef.te.lasso[i,,] <- coef.tr.lasso[which(tr.ind.o == te.ind.o[i]),,]
        }
        else if(belong1[i]){
          coef.te.lasso[i,,] <- coef.tr.lasso[which(tr.ind.o == (te.ind.o[i]-1)),,]
        }
        else if(te.ind.o[i] != 1) { # example: 5,6 in test set; 4, 7, 8, in training set
          coef.te.lasso[i,,] <- coef.te.lasso[i-1,,]
        }
      }
      
      ################### criteria 1: calculate the t(coef.diff) %*% solve(cov.te) %*% coef.diff for each lambda ######################
      ###############make a mistake here when calculating the criteria : modify it tomorrow
      Loss.te <- c()
      time.diff <- c(otpt.te[-1] - head(otpt.te, -1), 0)
      if(cv == 1){
        coef.te.whole <- sh.resp$coef[te.ind.o, ] #doesn't consider the intercept since the difference for this part is 0
        cov.te.solve <- lapply(1:n.te.o, function(i) solve(cov.te[-1,-1,i]))
        for(i in 1:length(Lambda)){
          coef.diff <- cbind(coef.te.lasso[,-1,i] - coef.te.whole[,-1]) # beta - \hat{beta}(from the whole sample)
          lsa <- lapply(1:n.te.o, function(j) t(coef.diff[j,]) %*% cov.te.solve[[j]] %*% coef.diff[j,] * time.diff[j])
          Loss.te <- c(Loss.te, sum(unlist(lsa))) # the loss vector for each subject at test data set
        }
      } else if(cv == 2){
        cov.score.te <- lapply(1:n.te.o, function(i) cov.score(sh.resp$coef, otpt.te[i], oz, oy, odelta, otpt)$cov)
        for(i in 1:length(Lambda)){
          score.list <- lapply(1:n.te.o, function(j) cov.score(coef.te.lasso[,,i], otpt.te[j], oz.te, oy.te, odel.te, otpt.te)$score)
          score <- lapply(1:n.te.o, function(j) colSums(score.list[[j]]))
          lsa <- lapply(1:n.te.o, function(j) t(score[[j]]) %*% solve(cov.score.te[[j]]) %*% score[[j]] * time.diff[j])
          Loss.te <- c(Loss.te, sum(unlist(lsa))) # the loss vector for each subject at test data set
        }
      }
      ##########the loss matrix at each lambda for each time point:difference form##########
      Loss.V <- rbind(Loss.V, Loss.te)
      #Loss.SD <- rbind(Loss.SD,apply(Loss.sd,2,sum)) 3used for 1se criteria
    }
  }
  if(count <= 2) {
    list(CV = 0)
  }
  else{
    Loss.sum <- apply(Loss.V,2,sum)/(count)
    Loss.sd <- apply(Loss.V,2,sd)/sqrt(count)
    min.ind <- which(Loss.sum == min(Loss.sum))
    cv.min <- Lambda[min.ind]
    se.ind <- which(Loss.sum >= (Loss.sum[min.ind]+Loss.sd[min.ind]))
    se.ind <- se.ind[which(se.ind >= min.ind)[1]] # to guarantee we search for a larger tunning parameter
    cv.1se <- Lambda[se.ind]
    list(loss = Loss.sum, cv.sd = Loss.sd, Loss.V=Loss.V, cv.min = cv.min, cv.1se = cv.1se, CV = 1)
  }
}


#################################################################################################
######################## score function at each time points for cv = 2###########################
cov.score <- function(coef, curr.tpt, z, y, d, otpt){
  n <- nrow(z); p <- ncol(z)
  ar.idx <- y >= curr.tpt # the index of subject at risk 
  fail.idx <- (y == curr.tpt) #the subject failed currently
  
  ### reconsruct the coeff with ties
  coef.whole <- c()
  for(i in 1:length(otpt)){
    ties <- sum(y == otpt[i]) #if ties = 1, it's the unique time point; otherwise, it have more than one
    if(ties == 1) coef.whole <- rbind(coef.whole, coef[i,])
    else coef.whole <- rbind(coef.whole, matrix(rep(coef[i,], ties), nrow = ties, byrow = T))
  }
  
  if(sum(fail.idx) == 1) coef.fail <- coef.whole[fail.idx,] 
  else coef.fail <- apply(coef.whole[fail.idx,],2,mean)
  #coef.fail <- ifelse(sum(fail.idx) == 1, coef.whole[fail.idx,], apply(coef.whole[fail.idx,],2,mean))
  
  counting.process <- (y <= curr.tpt & d == 1)# denote the value of N_i(t)
  
  
  if(sum(ar.idx)==n) martingale.process <-  as.vector(counting.process - exp(z %*% coef.fail))
  ###first failed subject
  else if(sum(ar.idx)==(n-1)){
    martingale.process.arrisk <-  counting.process[ar.idx] - exp(z[ar.idx,] %*% coef.fail)
    martingale.process.past <- counting.process[1-ar.idx] - exp(sum(z[1-ar.idx,] * coef.whole[1-ar.idx,]))
    martingale.process <- c(martingale.process.past,  martingale.process.arrisk)
  }
  else if(sum(ar.idx) > 1 & sum(ar.idx) < (n-1)){
    martingale.process.arrisk <-  counting.process[ar.idx] - exp(z[ar.idx,] %*% coef.fail)
    martingale.process.past <- counting.process[1-ar.idx] - exp(apply(z[1-ar.idx,] * coef.whole[1-ar.idx,], 1, sum))
    martingale.process <- c(martingale.process.past,  martingale.process.arrisk)
  }
  else {#only one subject left
    martingale.process.arrisk <-  counting.process[ar.idx] - exp(sum(z[ar.idx,] * coef.fail))
    martingale.process.past <- counting.process[1-ar.idx] - exp(apply(z[1-ar.idx,] * coef.whole[1-ar.idx,], 1, sum))
    martingale.process <- c(martingale.process.past,  martingale.process.arrisk)
  }
  
  
  score <- z * martingale.process/sqrt(n) ##### 1/root(n) * the summation of score
  score.cov <- t(score) %*% score ####the consistent est of cov of score
  list(cov=score.cov, score=score)
}

##########the value of least square approximation for score version ########## 
score.lsa <- function(cov, score){
  ic <- score %*% solve(cov) %*% t(score)
  sum(ic)
}