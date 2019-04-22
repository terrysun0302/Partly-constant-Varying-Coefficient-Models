#================================================================================================================#
#== this is an example of using the proposed method with an open source data set ================================#
#================================================================================================================#

#================================================================================================================#
#================== introduction motivation example: veteran data (open source) =================================#
#================================================================================================================#


#=============================================================#
#========== Step 1: import packages and R codes ==============#
#=============================================================#

### load important packages
library(genlasso)
library(survival)
library(dplyr)

### load all related R files
path <- paste0("C:/Users/hsun16/Desktop/")
dir <- paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/")
### Peng and Huang's method
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/path.R"))
### some important supporting functions for calculation in section 3.3
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/supportFunc.R"))
### cross validation function
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/cv.lsa.921.r"))
### resampling inference function
source(paste0(path, "/HSun/Research Project/Time varying model and variable selection/Fused lasso/final code/perturbed.inference.R"))

#=============================================================#
#========== Step 2: clean the data ===========================#
#=============================================================#

### clean the data
select <- dplyr::select
data(veteran)
veteran <- veteran %>% 
  select(celltype, time, status, karno) %>%
  mutate("sq" = ifelse(celltype == "squamous",1,0),"small" = ifelse(celltype == "smallcell",1,0), "adeno" = ifelse(celltype == "adeno",1,0))


veteran$celltype <-relevel(veteran$celltype, ref = "large")

### fit Cox model as benchmark and initial estimation for the proposed method.
mod <- coxph(Surv(time,status) ~ karno + celltype, data = veteran)

base <- log(basehaz(mod)[,1]); coeff <- coefficients(mod)
Initial <- cbind(base, matrix(rep(coeff,length(base)), length(base),length(coeff), byrow = T))
colnames(Initial) <- c("base","karno","sq","small","adeno")


data <- veteran %>% select(time, status, karno ,sq, small, adeno)
U <- data$time; Delta <- data$status; Z <- data %>% select(karno ,sq, small, adeno)


#==============================================================#
#====== Step 3: fit the proposed model ========================#
#==============================================================#
set.seed(02082018)

### fit the Peng and Huang's model (2007) and generate some important elements for later computation
vet1 <- cox_temporal_op(Z,U,Delta,init.fun=bini.fun, B=250, initial = Initial, loss = 2) #original fit
### fit the genlasso path using the Peng and Huang's estimation
whole.path.vet <- fit.genlasso(vet1,Design =T) #genlasso with Design = T

### cross validation to select the tunning parameter
set.seed(12345)
id1.vet <- cv.lsa2(vet1,CV=5,low=0.1, upp=0.9, V = F, cv = 1, Initial = Initial, loss = 2)
### the resampling method to do the covariane estimation
pi.vet <- perturbed.inference(vet1, id1.vet$cv.1se)

######### to find the indx of the coefficient in the whole.path
a <- match(vet1$time[11:75], vet1$time[vet1$unique.del==1])


#==============================================================#
#====== Step 4: generate the figure in section 1 ==============#
#==============================================================#
path <- paste0("C:/Users/hsun16/Desktop/HSun/Research Project/Time varying model and variable selection/Fused lasso")
pdf(paste0(path, "/final code/plot/introduction_data.pdf"), width = 10, height  = 10)

#par(mfrow = c(2,3))
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(3, 4, 2, 0), # two rows of text at the outer left and bottom margin
    mar = c(2, 1, 1, 1))
#######################################################################################################################
############################## the bias plot ##########################################################################
#######################################################################################################################
plot(vet1$time[11:75]/365, vet1$coef[11:75,2], type = 's', main = "(a)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.08, 0), col = 1)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[11:72], type = 's', col = 2)
legend("bottomleft", c( "PH(Card: 65)", "Proposed(Card: 14)"), lty = c(1,1), col = c(1,2))

plot(vet1$time[11:75]/365, vet1$coef[11:75,3], type = 's', main = "(b)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.5,1.2), col = 1)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97+11):(97+72)], type = 's', col = 2)
legend("bottomleft", c("PH(Card: 65)", "Proposed(Card: 10)"), lty = c(1,1), col = c(1,2))


plot(vet1$time[11:75]/365, vet1$coef[11:75,4], type = 's', main = "(c)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.5,2), col = 1)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97*2+11):(97*2+72)], type = 's', col = 2)
legend("bottomleft", c("PH(Card: 65)", "Proposed(Card: 5)"), lty = c(1,1), col = c(1,2))

plot(vet1$time[11:75]/365, vet1$coef[11:75,5], type = 's', main = "(d)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.5,2), col = 1)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97*3+11):(97*3+72)], type = 's', col = 2)
legend("bottomleft", c("PH(Card: 65)", "Proposed(Card: 5)"), lty = c(1,1), col = c(1,2))




#title(ylab = "Empirical Coverage Rate", outer = TRUE, line = 2, cex.lab = 2, adj = 0.1, font.lab = 2)
title(ylab = "Estimated Coefficient", outer = TRUE, line = 2, cex.lab = 2, adj = 0.5, font.lab = 2)
#title(ylab = "Empirical Bias", outer = TRUE, line = 2, cex.lab = 2, adj = 0.9,  font.lab = 2)
title(xlab = "Test time (years)", outer = TRUE, line = 2, cex.lab = 2)
#title(main = "Adverse events rate 45%", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.1)
#title(main = "Adverse events rate 20%", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.85)
#title(main = "Z2", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.85)
#title(main = expression(paste("Scenario 1, n = 150, ", beta[2], "(t) = -0.5")), outer = TRUE, line = 4.5, cex.main = 3.5, adj = 0.5)
##########################################################################################################################
############################################## end #######################################################################
##########################################################################################################################
dev.off()




### add the black_white version
path <- paste0("C:/Users/hsun16/Desktop/HSun/Research Project/Time varying model and variable selection/Fused lasso")
pdf(paste0(path, "/final code/plot/introduction_data_bw.pdf"), width = 10, height  = 10)

#par(mfrow = c(2,3))
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(3, 4, 2, 0), # two rows of text at the outer left and bottom margin
    mar = c(2, 1, 1, 1))
#######################################################################################################################
############################## the bias plot ##########################################################################
#######################################################################################################################
plot(vet1$time[11:75]/365, vet1$coef[11:75,2], type = 's', main = "(a)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.08, 0), col = "gray", lwd = 2)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[11:72], type = 's', col = 1, lwd = 2)
abline(a = mod$coefficients[1], b = 0, col = "gray", lty = 2, lwd = 2)
legend("bottomleft", c("Prop. Hazards" ,"Peng & Huang (2007)", "Proposed"), lty = c(2,1,1), col = c("gray","gray", 1), lwd = 2)

plot(vet1$time[11:75]/365, vet1$coef[11:75,3], type = 's', main = "(b)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.7,1.2), col = "gray", lwd = 2)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97+11):(97+72)], type = 's', col = 1, lwd = 2)
abline(a = mod$coefficients[2], b = 0, col = "gray", lty = 2, lwd = 2)
legend("bottomleft", c("Prop. Hazards" ,"Peng & Huang (2007)", "Proposed"), lty = c(2,1,1), col = c("gray","gray", 1), lwd = 2)


plot(vet1$time[11:75]/365, vet1$coef[11:75,4], type = 's', main = "(c)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.5,2), col = "gray", lwd = 2)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97*2+11):(97*2+72)], type = 's', col = 1, lwd = 2)
abline(a = mod$coefficients[3], b = 0, col = "gray", lty = 2, lwd = 2)
legend("bottomleft", c("Prop. Hazards" ,"Peng & Huang (2007)", "Proposed"), lty = c(2,1,1), col = c("gray","gray", 1), lwd = 2)

plot(vet1$time[11:75]/365, vet1$coef[11:75,5], type = 's', main = "(d)", xlab = "", ylab = "", frame.plot = F, ylim = c(-0.5,2), col = "gray", lwd = 2)
box(bty="l")
lines(vet1$time[11:75][!is.na(a)]/365, coef(whole.path.vet$out, lambda = id1.vet$cv.1se)$beta[(97*3+11):(97*3+72)], type = 's', col = 1, lwd = 2)
abline(a = mod$coefficients[4], b = 0, col = "gray", lty = 2, lwd = 2)
legend("bottomleft", c("Prop. Hazards" ,"Peng & Huang (2007)", "Proposed"), lty = c(2,1,1), col = c("gray","gray", 1), lwd = 2)



#title(ylab = "Empirical Coverage Rate", outer = TRUE, line = 2, cex.lab = 2, adj = 0.1, font.lab = 2)
title(ylab = "Estimated Coefficient", outer = TRUE, line = 2, cex.lab = 2, adj = 0.5, font.lab = 2)
#title(ylab = "Empirical Bias", outer = TRUE, line = 2, cex.lab = 2, adj = 0.9,  font.lab = 2)
title(xlab = "Test time (years)", outer = TRUE, line = 2, cex.lab = 2)
#title(main = "Adverse events rate 45%", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.1)
#title(main = "Adverse events rate 20%", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.85)
#title(main = "Z2", outer = TRUE, line = 1, cex.main = 2.5, adj = 0.85)
#title(main = expression(paste("Scenario 1, n = 150, ", beta[2], "(t) = -0.5")), outer = TRUE, line = 4.5, cex.main = 3.5, adj = 0.5)
##########################################################################################################################
############################################## end #######################################################################
##########################################################################################################################
dev.off()
