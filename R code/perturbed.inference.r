#=========================================================================#
#============= resampling method for covariance estimation ===============#
#=========================================================================#

library(Matrix)
library(genlasso)
path <- "C:/Users/hsun16/Desktop"
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/datgen.R"))
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/path.R"))
#source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/path.newton.R"))
#source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/path.newton2.R"))
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/path.loss2.R"))
#source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/cv.lsa.R"))
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/cv.lsa.921.R"))

###### set up the penalty matrix
###### set up the penalty matrix
perturbed.inference <- function(mod, lambda){
  if(!is.null(mod$sdmatrix)){
    B <- dim(mod$barray)[3]
    m.z <- apply(mod$oz,2,mean)
    n <- nrow(mod$oz);p <- ncol(mod$oz)
    risk <-  which(mod$unique.del==1)
    
    D = sd(mod$oz[,2]) * getD1d(length(risk))
    p.coef.mat <- array(numeric(0), dim = c(length(risk), p, B))
    
    if(p > 2){#have more than one covariate
      for(i in 3:p){
        D <- bdiag(D,sd(mod$oz[,i]) * getD1d(length(risk)))
      }
    }
    D <- as.matrix(D)
    
    Design.C <- getVY(mod$covmatrix[,,risk],length(risk),mod$coef[risk,-1], mod$time[risk])
    Y <- Design.C$Y
    X <- Design.C$Design
    out.gen <- genlasso(y=Y, X=X, D = D)###get the smoothed path for each lambda
    
    cov <- mod$covmatrix[,,mod$unique.del==1]
    coef.unpe <- mod$coef[mod$unique.del==1,]
    path.beta <- matrix(coef(out.gen, lambda = lambda)$beta, byrow = F, ncol = p-1)
    information <- lapply(1:sum(mod$unique.del), function(i) solve(cov[,,i])) 
    intercept.pe <- lapply(1:sum(mod$unique.del), function(i) coef.unpe[i,1] - 1/information[[i]][1,1] * t(information[[i]][1,-1]) %*% 
                             (path.beta[i,] - coef.unpe[i,-1]))
    intercept.pe <- unlist(intercept.pe) 
    
    beta.lambda.hat <- cbind(intercept.pe, path.beta)
    
    for(i in 1:B){
      Design.C <- getVY(mod$covmatrix[,,risk],length(risk),mod$barray[risk,-1,i], mod$time[risk])
      Y <- Design.C$Y
      X <- Design.C$Design
      out.gen <- genlasso(y=Y, X=X, D = D)###get the smoothed path for each lambda
      path.beta <- matrix(coef(out.gen, lambda = lambda)$beta, byrow = F, ncol = p-1)
      p.coef.mat[,-1,i] <- path.beta
      coef.unpe <- mod$barray[mod$unique.del==1,,i]
      intercept.pe <- lapply(1:sum(mod$unique.del), function(i) coef.unpe[i,1] - 1/information[[i]][1,1] * t(information[[i]][1,-1]) %*% 
                               (path.beta[i,] - coef.unpe[i,-1]))
      p.coef.mat[,1,i] <- unlist(intercept.pe)
    }
    sdmat <- apply(p.coef.mat, c(1,2), sd)
    ci.up <- beta.lambda.hat + 1.96 * sdmat
    ci.low <- beta.lambda.hat - 1.96 * sdmat
    ci.up.percen <- apply(p.coef.mat,c(1,2),function(i) quantile(i, prob = 0.975))
    ci.low.percen <- apply(p.coef.mat,c(1,2),function(i) quantile(i, prob = 0.025))
    
    ######### confidence band ##########
    #sweep:subtrat a matrix from an array
    J_star <- apply(abs(sweep(p.coef.mat,1:2,beta.lambda.hat))[-c(1:as.integer(n*0.1)),,],c(3,2), max) # doesn't the firs 10%
    r_alpha <- apply(J_star, 2, quantile, prob = 0.95)
    cbarray <- lapply(1:p, function(i) cbind(beta.lambda.hat[,i] - r_alpha[i], beta.lambda.hat[,i] + r_alpha[i]))
    
    list(sdmat = sdmat, ci.up.normal = ci.up, ci.low.normal = ci.low, ci.up.percen = ci.up.percen, ci.low.percen = ci.low.percen,
         cband = cbarray)
  }
  else return("unavailable")
}

# fp.10.2 <- list()
# count = 0
# for(i in 1:iter){
#   count <- count + 1
#   if(count <= 100) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[count]]
#   else if(count >= 101 & count < 200) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[count-101+2]]
#   else if(count >= 200 & count < 299) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[count-200+2]]
#   else if(count >= 299 & count < 398) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[count-299+2]]
#   else if(count >= 398 & count < 497) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[1]][[count-398+2]]
#   else if(count >= 497 & count < 596) sh.resp <- Result[[1]][[1]][[1]][[1]][[1]][[count-497+2]]
#   else if(count >= 596 & count < 695) sh.resp <- Result[[1]][[1]][[1]][[1]][[count-596+2]]
#   else if(count >= 695 & count < 794) sh.resp <- Result[[1]][[1]][[1]][[count-695+2]]
#   else if(count >= 794 & count < 893) sh.resp <- Result[[1]][[1]][[count-794+2]]
#   else if(count >= 893 & count < 992) sh.resp <- Result[[1]][[count-893+2]]
#   else if(count >= 992 & count < 1001) sh.resp <- Result[[count-992+2]]
#   
#   fp.10.2[[i]] <- sh.resp
#   #perturbed.ci[[i]] <- perturbed.inference(sh.resp, lambda = lambda.record[i,2])
# }
# save(fp.10.2, file = "result.rda")
# 
# 
# load("result.rda")
# load("lambda.record.rda")
# load("result.ws.rda")
# result <- fp.10.2
# iter <- 1000
# ##########Library for parallel computing##########
# #library(foreach)
# #library(doParallel)
# #no_cores <- detectCores() - 1
# #print(no_cores)
# #c1 <- makeCluster(no_cores)
# #registerDoParallel(c1)
# 
# # perturbed.ci <- list()
# # system.time(
# #   for(i in 1:iter){
# #     perturbed.ci[[i]] <- perturbed.inference(result[[i]], lambda.record[i,2])
# #   }
# # )
# 
# system.time(perturbed.ci <- foreach(i=1:iter, .combine = "list", .multicombine = T,.packages='genlasso') %dopar% {
#   #library(genlasso)
#   sh.resp <- result[[i]]
#   mod <- perturbed.inference(sh.resp, lambda.record[i,2])
#   mod
# } 
# )
# 
# perturbed.ci <- list()
# system.time(
#   for(i in 551:1000){
#     perturbed.ci[[i]] <- perturbed.inference(result[[i]], lambda.record[i,2])
#   }
# )
# 
# eta <- 0
# test.tp <- seq(0.15, 2.50, length.out = 50)##########the approximated 10% -- 90% event time
# true.0 <- log(test.tp) - log(exp(-3*test.tp+3)/(1+exp(-3*test.tp+3))+1.05)
# true.1 <- log(exp(-3*test.tp+3)/(1+exp(-3*test.tp+3))+1.05)
# true.2 <- rep(eta, length(test.tp))
# true.path <- cbind(true.1,true.2)
# 
# 
# setwd("C:/Users/hsun16/Desktop/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/output_ciband.250")
# # alloutputs <- list.files(path=".", pattern="job*")
# alloutputs <- paste0("job",1:1000,".rda")
# cover2.fuse1 <- array(numeric(0), dim = c(50, 4, 1000))
# sdmat2.fuse1 <- array(numeric(0), dim = c(50, 4, 1000))
# for(i in 1:1000){
#   load(alloutputs[i])
#   grid.idx <- apply(as.matrix(test.tp), 1, function(x){idx <- max(which(sh.resp$time[sh.resp$unique.del == 1] <= x)) })
#   fused.ci.normal.low <- perturbed.ci$ci.low.normal[grid.idx,]
#   fused.ci.normal.up <- perturbed.ci$ci.up.normal[grid.idx,]
#   fused.ci.percen.low <- perturbed.ci$ci.low.percen[grid.idx,]
#   fused.ci.percen.up <- perturbed.ci$ci.up.percen[grid.idx,]
#   sdmat.fuse1[,,i] <- perturbed.ci$sdmat[grid.idx,]
#   
#   cover2.fuse1[,1:2,i] <- true.path >= fused.ci.normal.low & true.path <= fused.ci.normal.up
#   cover2.fuse1[,3:4,i] <- true.path >= fused.ci.percen.low & true.path <= fused.ci.percen.up
# }
# cover2 <- apply(cover2.fuse1, c(1,2), mean)
# 
# 
# cover2.fuse <- array(numeric(0), dim = c(50, 4, 1000))
# sdmat2.fuse <- array(numeric(0), dim = c(50, 4, 1000))
# for(i in 1:1000){
#   grid.idx <- apply(as.matrix(test.tp), 1, function(x){idx <- max(which(result[[i]]$time[result[[i]]$unique.del == 1] <= x)) })
#   fused.ci.normal.low <- perturbed.ci[[i]]$ci.low.normal[grid.idx,]
#   fused.ci.normal.up <- perturbed.ci[[i]]$ci.up.normal[grid.idx,]
#   fused.ci.percen.low <- perturbed.ci[[i]]$ci.low.percen[grid.idx,]
#   fused.ci.percen.up <- perturbed.ci[[i]]$ci.up.percen[grid.idx,]
#   sdmat.fuse[,,i] <- perturbed.ci[[i]]$sdmat[grid.idx,]
#   
#   cover.fuse[,1:2,i] <- true.path >= fused.ci.normal.low & true.path <= fused.ci.normal.up
#   cover.fuse[,3:4,i] <- true.path >= fused.ci.percen.low & true.path <= fused.ci.percen.up
# }
# 
# cover <- apply(cover.fuse, c(1,2), mean)
# sdd.fuse <- apply(sdmat.fuse, c(1,2), mean)
