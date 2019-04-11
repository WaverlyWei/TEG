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
source((here("Lucy Project","R", "utils.R")))
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
complete_cov_1 <- read.csv("imputed_cov_11.csv")
complete_cov_1 <- complete_cov_1 %>% select(-X)

# impute on covariate set 2 
#imp_cov_2 <- amelia(cov_2)
#write.amelia(imp_cov_2, file.stem = "imputed_cov_2")
complete_cov_2 <- read.csv("imputed_cov_21.csv")
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
#outcome <- getMort6()

######################### Outcome 3 #####################
# death == 1 
# mortality 24h as outcome
#outcome <- getMort24()

######################### Outcome 4 #####################
# death == 1 
# mortality at discharge as outcome
outcome <- getMortDisch()
Y <- outcome[[1]]
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


