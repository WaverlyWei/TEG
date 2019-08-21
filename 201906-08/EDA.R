## Variable Importance Selection and CTMLE w/ hemostasis as outcome

rm(list=ls())
library(MASS)
library(SuperLearner)
library(origami)
library(glmnet)
library(randomForest)
library(ROCR)
library(cvAUC)
library(doParallel)
library(readstata13)
library(tmle)
library(xtable)
library(data.table)
library(R6)
library(tidyverse)
library(dplyr)
library(Amelia)
library(SIS)
library(ctmle)
library(vimp)
library(superheat)
library(pheatmap)
library(RColorBrewer)
library(heatmap.plus)
source((here("Desktop/Lucy Project","R", "utils.R")))
library(here)

set.seed(12313197)
options(digits = 4)
#setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Trauma and Coagulation (White Space Conflict)/ACIT/Data")
#setwd('/Users/waverlywei/Desktop/Lucy project')
dataf=read.dta13(here("Desktop/Lucy Project","ACIT_6Jan17_deidentified1.dta"))
plt0to6=as.numeric(dataf$icu_0to6h_plt_units>0)
bloodUnit0to6=as.numeric(dataf$icu_0to6h_blood_units>0)
ffp0to6=as.numeric(dataf$icu_0to6h_ffp_units>0)
dataf=data.frame(dataf,plt0to6,bloodUnit0to6,ffp0to6)
rm(plt0to6,bloodUnit0to6,ffp0to6)

# select covariates 
# split covariates into two sets
# make computation easier
cov_1 <- dataf %>% select(hr0_inr,hr0_ptt, sample0h_ex_ct, sample0h_crt_r, 
                        sample0h_ck_r, sample0h_crt_act, sample0h_ex_cft, 
                        sample0h_crt_k, sample0h_ck_k,sample0h_fib_mcf, 
                        sample0h_ff_ma, sample0h_ex_alpha, sample0h_crt_alpha, 
                        sample0h_ck_alpha)

# ffp, cryo: only have "icu_0to6h_units"
# plt: hr0_plts(predictor)  hr6_plts (outcome)
cov_2 <- dataf %>% select(sample0h_ex_a10, sample0h_ex_a20, sample0h_ex_mcf, 
                          sample0h_crt_ma, sample0h_ck_ma,sample0h_ex_ml, 
                          sample0h_ck_ly30, sample0h_ck_ly60,hr0_plts)
# char to numeric
cov_2$sample0h_ck_ly60 <- as.numeric(cov_2$sample0h_ck_ly60)

#note: use icu_0to6h_plts (outcome) 
out <- dataf %>% select(icu_0to6h_ffp_units,icu_0to6h_cryo_units,
                        icu_0to6h_plt_units,cyklokaproninpt)

out$cyklokaproninpt <- ifelse(out$cyklokaproninpt == "No",0,1)

# Imputation: prep for vimp 
# impute on covariate set 1
#imp_cov_1 <- amelia(cov_1)
#write.amelia(imp_cov_1, file.stem = "imputed_cov_1")
complete_cov_1 <- read.csv(here("Desktop/Lucy Project","imputed_cov_11.csv"))
complete_cov_1 <- complete_cov_1 %>% select(-X)

# impute on covariate set 2 
#imp_cov_2 <- amelia(cov_2)
#write.amelia(imp_cov_2, file.stem = "imputed_cov_2")
complete_cov_2 <- read.csv(here("Desktop/Lucy Project","imputed_cov_21.csv"))
complete_cov_2 <- complete_cov_2 %>% select(-X)

## -------------------------------------------------------------
## get variable importance!
## ------------------------------------------------------------

# remove NA outcomes 
# put this into loop 

#NA_index <- apply(df_all[,1:20],1, function(x) any(is.na(x)))
#df_all <- df_all[!NA_index,]

learner.lib <- c("SL.mean", "SL.xgboost", "SL.glmnet", "SL.randomForest")


# everytime create a new dataframe 
var_imp <- function(cov,df){
for (i in 1:4){
  for (j in 1:ncol(cov)){
    print(i)
    df_now <- df[,c(i,5:ncol(df))]
    
    # remove NA outcomes 
    NA_index <- apply(df_now[,1:ncol(df_now)],1, function(x) any(is.na(x)))
    df_now <- df_now[!NA_index,]
    
    co <- df_now[,2:ncol(df_now)]
    
    ## the full conditional mean
    full_regression <- SuperLearner(Y = df_now[,1], X = co, SL.library = learner.lib)
    full_fit <- full_regression$SL.predict
    ## the reduced conditional mean
    reduced_regression <- SuperLearner(Y = full_fit, X = co[,j, drop = FALSE], SL.library = learner.lib)
    reduced_fit <- reduced_regression$SL.predict
    
  
    ## get the variable importance estimate, SE, and CI
    vimp <- vimp_regression(Y = df_now[,1], f1 = full_fit, f2 = reduced_fit, indx = j, run_regression = FALSE)
    out_name <- c(out_name, names(df)[i])
    cov_name <- c(cov_name,names(co)[j])
    imp = c(imp,vimp$est)
    low = c(low, vimp$ci[1])
    high = c(high, vimp$ci[2])
  }
}
  return(list(imp,low,high,out_name,cov_name))
}


# run on covariate set 1
df_1 <-  cbind(out,complete_cov_1)
imp <- low <- high <- out_name <- cov_name <- NULL
result <- var_imp(complete_cov_1,df_1)

save(result,file = "vimp_1.Rdata")

# run on covariate set 2 
df_2 <- cbind(out,complete_cov_2)
imp <- low <- high <- out_name <- cov_name <- NULL

result.2 <- var_imp(complete_cov_2,df_2)

save(result.2,file = "vimp_2.Rdata")

## -------------------------------------------------------------
## Vimp Result Analysis
## ------------------------------------------------------------
load("vimp_1.Rdata")
df_vimp_1 <- data.frame(result[[4]],result[[5]],result[[1]], result[[2]],result[[3]])
names(df_vimp_1)<-c("out_name","cov_name","vimp","ci_low","ci_high")
  
# organize df into a matrix for plotting
temp.df.1 <- matrix(df_vimp_1$vimp, nrow = 4, ncol = 14,byrow = TRUE)
rownames(temp.df.1) <- names(out)
colnames(temp.df.1) <- names(complete_cov_1)


## Explore cov_set_2
load("vimp_2.Rdata")
df_vimp_2 <- data.frame(result.2[[4]],result.2[[5]],result.2[[1]], 
                         result.2[[2]],result.2[[3]])
names(df_vimp_2)<-c("out_name","cov_name","vimp","ci_low","ci_high")

# organize df into a matrix for plotting
temp.df.2 <- matrix(df_vimp_2$vimp, nrow = 4, ncol = 9,byrow = TRUE)
rownames(temp.df.2) <- names(out)
colnames(temp.df.2) <- names(complete_cov_2)

df_all <- cbind(temp.df.1,temp.df.2)

## Top 5 important variables for each blood product 
var_df <- matrix(NA,nrow = 4,ncol = 5)
for (i in 1:4){
  vars <- names(sort(df_all[i,],decreasing = TRUE)[1:5])
  var_df[i,] <- vars
}
row.names(var_df) <- names(out)
xtable(var_df)

# raw plot w/o clustering variable 

#-------------------------
#        Plot
#--------------------------
heatmap(df_all,cexRow = 1,
        main="Varibale Importance of Biomarkers VS. Blood Products",
        # draw vertical lines 
        add.expr = abline(v=c(3,5,9.5,12.5,16,19)))


pheatmap(df_all,
        main="Varibale Importance of Biomarkers VS. Blood Products")
        # draw vertical lines 
        #add.expr = abline(v=c(3,5,9.5,12.5,16,19)))
#------------------
# CTMLE
#-------------------

#outcome <- as.numeric(dataf$sample6h_ck_ly30)

######################### Outcome 1 #####################
## Hemostasis
# hemostatsis == 1 
#Y <- getY()


######################### Outcome 2 #####################
# death == 1 
# mortality 6h as outcome
outcome <- getMort6()

######################### Outcome 3 #####################
# death == 1 
# mortality 24h as outcome
#outcome <- getMort24()

######################### Outcome 4 #####################
# death == 1 
# mortality at discharge as outcome
#outcome <- getMortDisch()
Y <- outcome[[1]]
#Y <- as.character(unlist(Y))
#Y <- ifelse(Y == "Dead",1,0)
names(Y) <- "Y"
index <- outcome[[2]]

impute_W <- getW()
#baseline <- getY(complete_cov_2$sample0h_ck_ly30,0,8)
#baseline <- complete_cov_1$sample0h_crt_k
#final_W <- cbind(impute_W,baseline)
final_W <- impute_W[index,]
#A <- getA("crt_alpha","icu_0to6h_cryo_units")[index]
A <- getA("crt_ma","icu_0to6h_plt_units")[index]
#A <- getA("crt_act","icu_0to6h_ffp_units")[index]

whole <- cbind(Y,final_W,A)





# whole <- whole %>% filter(!is.na(Y) &
#                             !is.na(A) & !is.na(final_W)) 

whole <- whole %>% filter(!is.na(Y) &
                            !is.na(A)) 



# w/o prespecifying Q  
ctmle_discrete_fit2 <- ctmleDiscrete(Y = whole$Y, A = whole$A, W = whole[,3:18] ,
                                     preOrder = FALSE, detailed = TRUE)
ctmle_discrete_fit2$est
ctmle_discrete_fit2$CI
ctmle_discrete_fit2$pvalue

# ======================
# Unadjusted ATE
# ======================
#unadjusted ate
t.1 = whole[whole$A==1,]
t.0 = whole[whole$A==0,]
# on/hemo == 1 
on_hemo <- sum(t.1$Y)/nrow(t.1)
# off/hemo == 1 
off_hemo <- sum(t.0$Y/nrow(t.0))
# ate
on_hemo - off_hemo

## Bar Plot 
# hemostasis as outcome
df1 <- data.frame(Treatment = factor(c("plt","plasma","cryo")),
                 ATE = c(0.1328,0.1111,0.2026),
                 lower = c(0.04296,0.09142,0.1042),
                 upper = c(0.22260,0.13073,0.3009))
p1 <- ggplot(df1, aes(Treatment, ATE, colour = Treatment)) +
     geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.8,fatten = 0.8) +
     theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
       ggtitle("ATE in Hemostasis Evaluated with CTMLE")+
  theme(plot.title = element_text(hjust = 0.5,size = 30),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  ylim(-0.5,0.5) + 
  geom_hline(yintercept = 0) 
p1

# mortality at 6h 
df2 <- data.frame(Treatment = factor(c("plt","plasma","cryo")),
                  ATE = c(-0.1365,-0.1692,-0.1584),
                  lower = c(-0.3239,-0.2159,-0.23466),
                  upper = c(0.0508,-0.1224,-0.08216))
p2 <- ggplot(df2, aes(Treatment, ATE, colour = Treatment)) +
  geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("ATE in Mortality at 6h Evaluated with CTMLE")+
  theme(plot.title = element_text(hjust = 0.5,size = 30),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  ylim(-0.5,0.5) + 
  geom_hline(yintercept = 0)
p2

## P-value plot 
# hemostasis
df3 <- data.frame(Treatment = factor(c("plt","plasma","cryo")),
                  pDR = c(0.271,0.2115,0.314),
                  pCTMLE = c(0.003764,0.00001673,0.00005436))


p3 <- ggplot(df3, aes(Treatment, pDR,colour = "DR",group = "pDR")) + 
       geom_line(size = 2) + geom_point(size = 3) +
       geom_line(aes(Treatment, pCTMLE,colour = "CTMLE",group = "pCTMLE"),size = 2)+
       geom_point(aes(Treatment, pCTMLE,colour = "CTMLE",group = "pCTMLE"),size = 3) + 
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("P-values with Doubly Robust Estimator vs.CTMLE (Hemostasis)")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) + 
  ylim(0,0.4) + ylab("P-value")

p3


# p-value : mortality as outcome 
df4 <- data.frame(Treatment = factor(c("plt","plasma","cryo")),
                  pDR = c(0.2598,0.1672,0.1531),
                  pCTMLE = c(0.01531,0.00001361,0.00004665))
p4 <- ggplot(df4, aes(Treatment, pDR,colour = "DR",group = "pDR")) + 
  geom_line(size = 2) + geom_point(size = 3) +
  geom_line(aes(Treatment, pCTMLE,colour = "CTMLE",group = "pCTMLE"),size = 2)+
  geom_point(aes(Treatment, pCTMLE,colour = "CTMLE",group = "pCTMLE"),size = 3) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("P-values with Doubly Robust Estimator vs.CTMLE (Mortality 6h)")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) + 
  ylim(0,0.4) + ylab("P-value")

p4

## Plot sparse data
cryo_dat <- dataf$icu_0to6h_cryo_units
plt_dat <- dataf$icu_0to6h_plt_units
plasma_dat <- dataf$icu_0to6h_ffp_units

sparse_dat_w <- data.frame(cryo_dat, plt_dat, plasma_dat)

sparse_dat <- melt(sparse_dat)
names(sparse_dat) <- c("Treatment", "Value")



library(naniar)
# using  geom_miss_point()
# sparse_p <- ggplot(sparse_dat,
#        aes(x = Treatment,
#            y = Value)) +
#   geom_miss_point(alpha = 0.1,size = 5) +
#   theme(panel.grid.major = element_blank(), 
#                                panel.grid.minor = element_blank(),
#                                panel.background = element_blank(), 
#                                axis.line = element_line(colour = "black")) 



sparse_p <- sparse_dat %>%
  bind_shadow() %>%
  ggplot(aes(x = Treatment,
             fill = Value_NA)) + 
  geom_density(alpha = 0.5) +  
  theme(panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     panel.background = element_blank(), 
                                     axis.line = element_line(colour = "black")) 


# simulation table
bias <- c(2766.8,10.8,1.3,0.4,0.3)
se <- c(22.61,13.52,11.05,10.41,10.41)
mse <- c(7706.3,18.4,12.2,10.8,10.6)
df <- cbind(bias,se,mse)
rownames(df) <- c("unadj","AIPTW","TMLE","greedy-CTMLE","SL-CTMLE")
colnames(df) <- c("bias(10E-3)","se(10E-2)","MSE(10E-3)")
library(xtable)
xtable(df)


## Simulation 
# Highly correlated A & W's 
library(locfit)
n <- 1000
W1 <- rbinom(n,1,0.5)
W2 <- rbinom(n,1,0.3)
W3 <- rbinom(n,1,0.01)
W4 <- rbinom(n,1,0.2)
W5 <- rbinom(n,1,0.8)
W6 <- rbinom(n,1,0.1)

W <- data.frame(W1,W2,W3,W4,W5,W6)

g0 <- expit(0.1*W1 - 0.3*W2 + 0.5*W3 + 0.05 * W4 + 0.2*W5 - 0.1*W6)
A <- ifelse(g0 >= 0.5, 1,0)

Y <- 5 + A + W1 + 3*W2 + 4*W3 + W4 + W5 + W6
Y <- ifelse(Y>=mean(Y),1,0)

dat <- data.frame(W1,W2,W3,W4,W5,W6,A,Y)

Q0 <- mean(dat[A==1,]$Y) - mean(dat[A==0,]$Y)

learner.lib <- c("SL.mean", "SL.xgboost", "SL.glmnet", "SL.randomForest")

# 100 simulations
tmle_res <- tmle(Y = Y, A = A,
                 W = W)
#tmle_res$estimates$ATE$psi
tmle_est <- tmle_res$estimates$ATE$psi
tmle_lower <- tmle_res$estimates$ATE$CI[1]
tmle_upper <- tmle_res$estimates$ATE$CI[2]

ctmle_res <- ctmleDiscrete(Y = Y, A = A, W = W ,
                           preOrder = TRUE, detailed = TRUE)
ctmle_est <- ctmle_res$est
ctmle_lower <- ctmle_res$CI[1]
ctmle_upper <- ctmle_res$CI[2]



# Highly correlated covariates 
df1 <- data.frame(Estimator = factor(c("TMLE","CTMLE")),
                  ATE = c(tmle_est,-ctmle_est),
                  lower = c(tmle_upper,-ctmle_upper),
                  upper = c(tmle_lower,-ctmle_lower))
p1 <- ggplot(df1, aes(Estimator, ATE, colour = Estimator)) +
  geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.8,fatten = 0.8) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Estimator Simulation Performance I")+
  theme(plot.title = element_text(hjust = 0.5,size = 30),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+ ylim(-1,0) +
 geom_hline(yintercept = Q0) 

p1






