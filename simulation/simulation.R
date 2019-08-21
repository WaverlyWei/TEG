library(here)
library(locfit)
library(tmle)
library(SuperLearner)
library(ctmle)
library(ggplot2)
library(microbenchmark)

set.seed(39172681)

# simulate high dimensional W's and highly correlated W&A
# 20 W's for now 
nCov <- 20
nObs <- 10
prob <- runif(nCov, 0.1,0.9)
Ws <- matrix(NA,nrow = nObs,ncol = nCov)

set.seed(39172681)
for (i in 1:nCov){
  W <- rbinom(nObs,1,prob[i])
  Ws[,i] <- W
}

A <- expit(0.5 * Ws[,3] + 0.19*Ws[,5] +
             0.1*Ws[,12] + 0.4*Ws[,20] -
             0.22*Ws[,18] - 0.52*Ws[,4] -
             0.37 * Ws[,7] - 0.19 * Ws[,19])
A <- ifelse(A>=0.5,1,0)
Y <- abs(A * Ws[,1] * Ws[,2] + Ws[,3] * Ws[,18] - Ws[,4] -
            Ws[,5]  + Ws[,6])

#Y <- rbinom(10,1,0.5)
varName <- c(paste("W", 1:nCov,sep = ""),"A","Y")

df <- data.frame(Ws,A,Y)
names(df) <- varName

true_ATE <- mean(df[A==1,]$Y) - mean(df[A==0,]$Y)
true_ATE
#[1] 0.125

# ======= parametric bootstrap =========== #
# theta* is ATE's variation
SL.lib <- c("SL.glm","SL.ridge","SL.mean","SL.xgboost")

# TMLE
tmle_est <- tmle(Y = df$Y, A = df$A, W = df[,1:20],
                 Q.SL.library = SL.lib,g.SL.library = SL.lib)
tmle_ate <- tmle_est$estimates$ATE$psi
tmle_sd <- tmle_est$estimates$ATE$var.psi

# CTMLE
# Bug: 
# 1. not because Y is binary, rbinom actually works 
set.seed(39172681)
ctmle_est <- ctmleDiscrete(Y = Y, A = df$A,W = df[,1:20],
                           #SL.library = SL.lib,
                           preOrder = TRUE, detailed = TRUE)
ctmle_ate <- ctmle_est$est
ctmle_sd <- ctmle_est$var.psi






# AIPW
calc_aipw <- function(d, Qform, g1W) {
   Q <- glm(Qform, data = d)
   QAW.pred <- predict(Q)
   # delete original A at column 21
   Q1W.pred <- predict(Q, newdata = data.frame(d[,-21], A = 1))
   Q0W.pred <- predict(Q, newdata = data.frame(d[,-21], A = 0))
   h <- d$A/g1W - (1- d$A)/(1 - g1W)
   return(list(psi = mean(h*(d$Y - QAW.pred) + Q1W.pred - Q0W.pred), 
          var = var(h*(d$Y - QAW.pred) + Q1W.pred - Q0W.pred)))
}



g <- glm(A ~ ., data = df[,-ncol(df)], family = "binomial")
g1W <- predict(g, type = "response")
wt <- A/g1W + (1 - A)/(1 - g1W)
wt.stab <- (A * mean(A) + (1 - A) * (1 - mean(A))) * wt
# IPW MSM
ipw.msm <- coef(glm(Y ~ A, data = df,weights = wt))
ipw_ate <- ipw.msm[2]
# IPW STAB.MSM
ipw.stab.msm <- coef(glm(Y ~ A, data = df,weights = wt.stab))
ipw.stab_ate <- ipw.stab.msm[2]
# aipw 
aipw <- calc_aipw(df, Qform = "Y ~ A", g1W = g1W)
aipw_ate <- aipw$psi
aipw_sd <- aipw$var


    
ParBootstrap <- function(ate,sd){
  theta_star <- vector()
  for (i in 1:1000){
    theta_star[i] <- mean(rnorm(length(df),
                                mean = ate,sd = sqrt(sd)))
  }
  theta_boot <- mean(theta_star) # bootstrap estimate of theta_hat
  boot_se <- sd(theta_star) # standard eorrs of the estimate
  bias <- theta_boot - true_ATE
  MSE <- mean((theta_boot- true_ATE)^2)
  CI <-c(theta_boot-1.96*boot_se,theta_boot +1.96*boot_se)
  return(list(mean = theta_boot,se = boot_se,bias = bias,mse = MSE,ci = CI))
}

# TO DO: 
# 1. Record simulation results 
# 2. Plot 
sim_res <- matrix(NA,nrow = 5, ncol = 6)
aipw_boot <- ParBootstrap(aipw_ate,aipw_sd)
sim_res[1,] <- unlist(aipw_boot)
tmle_boot <- ParBootstrap(tmle_ate,tmle_sd)
sim_res[2,] <- unlist(tmle_boot)
ctmle_boot <- ParBootstrap(ctmle_ate,ctmle_sd)
sim_res[3,] <- unlist(ctmle_boot)
ipw_boot <- ParBootstrap(ipw_ate,ctmle_sd)
sim_res[4,] <- unlist(ipw_boot)
ipw.stab_boot <- ParBootstrap(ipw.stab_ate,ctmle_sd)
sim_res[5,] <- unlist(ipw.stab_boot)

sim_res <- data.frame(sim_res)
names(sim_res) <- c("ATE","se","bias","mse","cilow","cihi")
rownames(sim_res) <- c("AIPW","IPW","SL-CTMLE","TMLE","Stablized IPW")
sim_res$estimator <- factor(rownames(sim_res))

# plot
pdf("/Users/waverlywei/Desktop/Thesis/simulation2.pdf") 
pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(sim_res, aes(x=estimator, y=ATE, group = 1,colour=estimator)) + 
  geom_errorbar(aes(ymin=cilow, ymax=cihi), width=.2, position=pd) +
  geom_point(position=pd,size = 3) + 
  ggtitle("Comparison of Estimator Performance")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0.01, linetype="dashed", color = "black")

dev.off()

xtable(sim_res[,-7])

# Time complexity: microbenchmark 
# example
N <- 100
p <- 100
Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]
tau <- 2
gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
W <- as.matrix(Wmat)
g <- 1/(1+exp(W%*%gcoef /3))
A <- rbinom(N, 1, prob = g)
epsilon <-rnorm(N, 0, 1)
Y  <- beta0 + tau * A + epsilon
# With initial estimate of Q
Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))
microbenchmark(ctmleDiscrete(Y = Y,
                                     A = A,
                                     W = W,
                                     Q = Q,
                                    preOrder = FALSE,
                                     detailed = TRUE),times = 1)

