#01/02/2017 re-organize the cox_temporal code for the resampling procedure proposed in Peng(2007)
#Use optim to solve the first failure time point(explain the reason to Brent)
#test the results using optim and newton-raphson, respectively.(show the table to Brent)

### 04/10/2019
#******************************************************************************#
# path function                                                      #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  z            design matrix                                                  #
#                                                                              #
#                                                                              #
#  y            the observed survival time                                     #
#                                                                              #
#                                                                              #
#  delta        censoring indicator                                            #
#                                                                              #
#  init.fun     initialization function for parameter estimation. Used in      #
#               simulation study. In real data example, we use 'initial'       #
#               argument                                                       #
#                                                                              #
#  B            the number of resampling for covariance estimation             #
#                                                                              #
#  initial      initialization parameter estimation. Usually we use the        #
#               parameter estimation of cox model here. See the veteran_cox.r  #
#               as an example.                                                 #
#  loss           the method to test DEGs. "1" is Gaussian mixture model; "2"  #
#                 is Anderson-darling normal test. The default value is "1"    #
#                                                                              #
#******************************************************************************#
cox_temporal_op <- function(z,y,delta,init.fun,B=0, initial = NULL, loss = 1){
  ###############################The optim function for the first failure time point##########
  fn1 <- function(beta,curr.tpt,z,y,d,og, loss) {
    #cat("curr.I = ",curr.I,"\n")
    ar.idx <- y >= curr.tpt # the index of subject at risk 
    fail.idx <- (y==curr.tpt & d==1) # the index of failure subject
    Nj <- sum(fail.idx) # the number of failure at a failure time(due with ties)
    zar <- z[ar.idx,] # z at-risk, say nar-by-p
    zarb <- as.vector(zar %*% beta) # nar-by-1
    ezarb <- exp(zarb) # nar-by-1
    
    if(loss == 1){
      if(Nj<2) S1 <- (1 + og[fail.idx])* z[fail.idx,]
      else S1 <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
      f1 <- (sum((S1 - t(zar) %*% ezarb)^2))
    }
    else if(loss == 2){
      if(Nj<2) fail.z <- (1 + og[fail.idx]) * z[fail.idx,]
      else fail.z <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
      f1 <- sum(ezarb) - t(beta) %*% fail.z
    }
    f1
  }
  
  # fn <- function(beta,curr.tpt,z,y,d,og) {
  #   #cat("curr.I = ",curr.I,"\n")
  #   ar.idx <- y >= curr.tpt # the index of subject at risk
  #   fail.idx <- (y==curr.tpt & d==1) # the index of failure subject
  #   Nj <- sum(fail.idx) # the number of failure at a failure time(due with ties)
  #   zar <- z[ar.idx,] # z at-risk, say nar-by-p
  #   zarb <- as.vector(zar %*% beta) # nar-by-1
  #   ezarb <- exp(zarb) # nar-by-1
  #   if(Nj<2) S1 <- (1 + og[fail.idx])* z[fail.idx,]
  #   else S1 <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
  #   
  #   f1 <- sum((S1 - t(zar) %*% ezarb)^2)
  #   
  #   if(Nj<2) fail.z <- (1 + og[fail.idx]) * z[fail.idx,]
  #   else fail.z <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
  #   f2 <- (sum(ezarb) - t(beta) %*% fail.z)
  #   
  #   #f1 <- sum(ezarb) - t(beta) %*% z[fail.idx,]
  #   s1 <- S1 - t(zar) %*% ezarb
  #   list(f1, f2, s1)
  # }
  
  BetaFITTER.first.optim <- function(binit,tpt,oz,oy,od,og=rep(0,nrow(z)), loss) { 
    # assume z is ordered, with intercept in 1st col.
    # objective function for beta(x[1])
    ft <- optim(f=fn1,p=binit,curr.tpt=tpt,z=oz,y=oy,d=od,og=og,loss=loss,control=list(trace=FALSE,maxit = 2000))
    list(b=ft$par,conv=ft$convergence)#conv = 0 means convergent
  }
  ##########################################################################################
  
  ###################################optim method for the temporal coefficient##############
  # objective function for beta(t), t>1 
  fn2 <- function(beta,curr.tpt,z,y,d,beta.last,og, loss) {
    ar.idx <- y >= curr.tpt # the index of subject at risk 
    fail.idx <- (y==curr.tpt & d==1) # the index of failure subject
    Nj <- sum(fail.idx) # the number of failure at a failure time(due with ties)
    zar <- z[ar.idx,] # z at-risk, say nar-by-p
    zarb <- as.vector(zar %*% beta) # nar-by-1
    ezarb <- exp(zarb) # nar-by-1
    ezarpb <- exp(as.vector(zar %*% beta.last))
    
    if(loss == 1){
      wt <- ezarb - ezarpb
      if(Nj<2) S1 <- (1 + og[fail.idx])* z[fail.idx,]
      else S1 <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
  
      if(sum(ar.idx)!=1) f1 <- (sum((-S1 + t(zar) %*% wt)^2))
      else f1 <- (sum((-S1 + t(zar) * wt)^2))
    }
    else if(loss == 2){
      betaz <- zar %*% beta
      if(Nj<2) fail.z <-  (1 + og[fail.idx]) * z[fail.idx,]
      else fail.z <- apply((1 + og[fail.idx]) * z[fail.idx,],2,sum)/Nj
       f1 <- sum(ezarb - ezarpb * betaz) - t(beta) %*% fail.z
    }
    
    f1
  }
  
  BetaFITTER.optim.sh <- function(binit,tpt,oz,oy,od,beta.last,og=rep(0,nrow(z)),loss) {
    
    ft <- optim(f=fn2,p=binit,curr.tpt=tpt,z=oz,y=oy,d=od,beta.last=beta.last,og=og,loss=loss,
                control=list(trace=FALSE,maxit = 2000))
    
    list(b=ft$par,conv=ft$convergence)#1, the solution; !1, some problems
  }
  #########################################################################################
  
  
  ###########################This function to estimate the time varying coeff##############
  PathFITTER.sh <- function(oz,oy,odelta,otpt,og=rep(0,nrow(z)),loss) {
    nbeta <- length(otpt) # the number of distinct time points
    ndata <- dim(oz)[1]#the number of subjects
    p <- dim(oz)[2] # the number of covariates
    bmat <- matrix(numeric(0),nbeta,p)#save the beta matrix for the estimation
    beta.next <- c(log(1/nbeta),rep(0,p-1)) # initial beta
    cnvrg.yesno <- rep(0,nbeta) # to record the convergence status
    del.yesno <- rep(0, nbeta)
    
    cntr <- 0
    p.tpt <- 0
    conv <- I <- 0
    while((I < nbeta) & (conv == 0)){
      #calcullate the failure index
      I <- I +1
      curr.tpt <- otpt[I] # we only consider unique time points in this path algorithm
      ar.idx <- oy >= curr.tpt
      fail.idx <- (oy==curr.tpt & odelta==1)
      Nj <- sum(fail.idx)
      #cat(I,": b[",oy[I],"] = ",round(beta.next,3),"\n")
      beta.last <- beta.next
      converge <- 0
      
      if ((Nj < 0.5)|(curr.tpt==p.tpt)){# if censored or the ties time points
        del.yesno[I] <- 0
        beta.next <- beta.last
        converge <- 2 				# continue due to censoring
      }
      else { 	# if failure
        #########Question : how to select the initial value
        del.yesno[I] <- 1
        #binit <- init.fun(curr.tpt)
        #binit <- rep(0,p)
        #binit <- c(-5, 0.5)
        binit <- initial[I,]
        #binit <- beta.last
        cntr <- cntr + 1 #count the number of failure time
        if (cntr==1) { 				    # first failure(important part)
          ft <- BetaFITTER.first.optim(binit,curr.tpt,oz,oy,odelta,og,loss)
          converge <- ft$conv
          if(converge == 0) beta.next <- ft$b
          else {
            conv <- 1
            break
          }
        }
        else { # cntr > 1 not the first failure
          if ((nbeta-I)>=0) { # run all the way to the end(until at least 1 subject left.)
            ft <- BetaFITTER.optim.sh(binit,curr.tpt,oz,oy,odelta,beta.last,og,loss)
            converge <- ft$conv #0 or 5 ok; 1, stop
            if(converge == 0) beta.next <- ft$b
            else {
              conv <- 1
              break
            }
          }
          else {     # no subject left and save it as the last eatimate
            beta.next <- beta.last 
            converge <- 3 
          }
          
        } # end loop for cntr >1
      } # end loop for failures
      bmat[I,] <- beta.next
      cnvrg.yesno[I] <- converge
      p.tpt <- curr.tpt
    } # end for 'I' loop
    list(path=bmat,converge=cnvrg.yesno, conv=conv, unique.del = del.yesno)
  }
  ####################################################################################
  ### End of definitions #############################################################
  ####################################################################################
  
  
  ####################################################################################
  ### Start computation  #############################################################
  ####################################################################################
  
  # append intercept, if necessary
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(is.matrix(z)&&(!any(z[,1]!=1))) z <- z
  else z <- cbind(rep(1,nrow(z)),as.matrix(z))
  
  dimz <- dim(z)
  #n <- dimz[1]
  #p <- dimz[2]
  
  # sort data
  # =-=-=-=-=-=-=-=
  order.time <- order(y)
  oy <- y[order.time]
  odelta <- delta[order.time]
  otpt <- unique(oy)
  #f.tpt <- 
  oz <- z[order.time,]
  
  ndata <- dimz[1]; nbeta <- length(otpt)
  p <- dimz[2]
  
  # the beta estimation for the original data
  ft0 <- PathFITTER.sh(oz,oy,odelta,otpt,og=rep(0,ndata),loss = loss)
  coef <- ft0$path
  unique.del <- ft0$unique.del
  
  m.z <- apply(oz, 2, mean)
  coef.center <- coef
  coef.center[,1] <- coef %*%  m.z

  # initialize, matrix/vector storage, when B > 0, we save the resampling results here.
  barray <- ciarray <- cbarray <- gmat <- sdarray <- st.coef <- st.sd <- sdmatrix <- covmatrix <- NULL
  count <- 0; betamatrix <- NULL
  
  if((B>0) & (ft0$conv == 0)) { #means this simulated dataset is solvable 
    #the est array for B resampling times
    barray <- barray.center <- array(numeric(0),dim=c(nbeta,p,B))
    ciarray <- cbarray <- array(numeric(0),dim=c(nbeta,p,2)) #confidence interval array
    #bmnext <- matrix(rep(beta.next,B),p,B,byrow=FALSE)
    
    #generate the beta array for B resampling times
    b <- 0
    while(b <= B){
      count <- count + 1
      ft_resp <- PathFITTER.sh(oz,oy,odelta,otpt,og=rnorm(ndata),loss=loss)
      if(ft_resp$conv == 0) {
        barray[,,b] <- ft_resp$path
        barray.center[,,b] <- barray[,,b]
        barray.center[,1,b] <- ft_resp$path %*% m.z
        b <- b + 1
      }
    }

    ##########################construct the covariance structure for beta(t)(pn*pn)#######################
    sdmatrix <- apply(barray, c(1,2), sd) # n-by-p matrix
    ######rewrite the resampling estimate which is easier to compute the covariance matrix
    
    ###### 10/16/2017: modify it here because I want to calculate the centralized covariance matrix
    betamatrix <- matrix(numeric(0), p*nbeta, B)
    index <- 1
    for(i in 1:nbeta){
      for(j in 1:p){
        betamatrix[index,] <- barray.center[i,j, ]
        index = index + 1
      }
    }
    ######calculate the centerlized-coefficient covariance matrix of coefficients at each time point 
    covmatrix <- array(numeric(0),dim = c(p,p,nbeta))
    for(i in 1:nbeta){
      covmatrix[,,i] <- cov(t(betamatrix[(i*p-p+1):(i*p),]))
    }
    
    #######################construct the confidence band for beta(t)#######################
    #sweep:subtrat a matrix from an array
    # J_star <- apply(abs(sweep(barray,1:2,sh.resp$coef))[-c(1:as.integer(n*0.2)),,],c(3,2), max) #say, B-by-p
    # r_alpha <- apply(J_star, 2, quantile, prob = 0.95)
    # cbarray[,,1] <- sweep(ft0$path,2,r_alpha)
    # cbarray[,,2] <- sweep(ft0$path,2,-r_alpha)
  }
  
  
  list(time=otpt,delta=odelta,unique.del=unique.del,oz=oz,oy=oy,converge=ft0$converge,coef=ft0$path,coef.center = coef.center ,sdmatrix=sdmatrix,
       barray=barray, count=count, conv=ft0$conv, covmatrix=covmatrix, betamatrix=betamatrix)
}