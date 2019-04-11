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

set.seed(12313197)
options(digits = 4)
#setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Trauma and Coagulation (White Space Conflict)/ACIT/Data")
setwd('/Users/weilinqing/Desktop/research project')
dataf=read.dta13("ACIT_6Jan17_deidentified1.dta")
plt0to6=as.numeric(dataf$icu_0to6h_plt_units>0)
bloodUnit0to6=as.numeric(dataf$icu_0to6h_blood_units>0)
ffp0to6=as.numeric(dataf$icu_0to6h_ffp_units>0)
dataf=data.frame(dataf,plt0to6,bloodUnit0to6,ffp0to6)
rm(plt0to6,bloodUnit0to6,ffp0to6)

# locate covariates 
nm = names(dataf)
# check position 
nm[grep("iss",nm)]
grep("hr0_basedefexc",nm)
grep("male",nm)
grep("age",nm)
grep("mechanism",nm)
grep("sbp",nm)
nm[grep("sbp",nm)]
grep("hr0_hr",nm)
grep("race",nm)
# iss
iss = dataf[,516]
# hr0_basexc
basexc = dataf[,26]
# gender
gender = as.character(dataf[,490])
gender[gender=="Male"] = 1
gender[gender=="Female"] = 0
gender = as.numeric(gender)
# age
age = dataf[,491]
# Mechanism ??
# sbp
sbp = dataf[,20]
# Hr
Hr = dataf[,18]
# race
race = as.character(dataf[,493])
unique(race)
race[race=="Unknown"] = 0
race[race=="White"] = 1
race[race=="Black"] = 2
race[race=="Asian"] = 3
race[race=="Native American"] = 4
race[race=="Pacific Islander"] = 5
race[race=="Other"] = 6
race = as.numeric(race)
W = cbind(iss,basexc,gender,age,sbp,Hr,race)


# clean up outcome vars
outcome = dataf %>% dplyr::select(contains("mortality"))
# death =1  alive = 0, for now, not remove the death ones
out = data.frame(outcome[,c(4,3)])
#colnames(out) = "mortality6h"
out$id = 1:nrow(out)
alive = out[out$mortalityat6h==0,]$id

## Targeting on ACT
act = dataf %>% dplyr::select(contains("crt_act"))
plasma = dataf %>% dplyr::select(contains("ffp_units"))

# remove blank labs
# nb = no blank labs at 2h, 4h 
nb = !is.na(act[alive,]$sample2h_crt_act) & !is.na(act[alive,]$sample4h_crt_act)
# valid row number
valid_obs = out[nb,]$id
head(valid_obs)

# clean-up everything
# NOTE: use mortality24h, otherwise no death 
final_out = out[valid_obs,]$mortalityat24h
final_act = act[valid_obs,1]
final_W = W[valid_obs,]
final_plasma = plasma[valid_obs,1]

# code treatment A
whole = cbind(final_out,final_act,final_plasma,final_W)
# comply
comply = whole[(whole[,2] >= 128 & !is.na(whole[,2])),]
nrow(comply)

# comply + plasma
get_plasma = comply[(comply[,3]>0),]
nrow(get_plasma)
# comply + no plasma
no_plasma = comply[(comply[,3]==0),]
nrow(no_plasma)

# not comply
not_comply = whole[(whole[,2] < 128 | is.na(whole[,2])),] 
nrow(not_comply)
# not comply + plasma
get_plasma_nc = not_comply[(not_comply[,3]>0),]
nrow(get_plasma_nc)
# not comply + no plt
no_plasma_nc = not_comply[(not_comply[,3]==0),]
nrow(no_plasma_nc)

# death rate?? maybe not remove deaths?? CORRECT THIS 

sum(get_plasma[,1],na.rm = TRUE) / nrow(get_plasma)
sum(no_plasma[,1],na.rm = TRUE) / nrow(no_plasma)
sum(get_plasma_nc[,1],na.rm = TRUE) / nrow(get_plasma_nc)
sum(no_plasma_nc[,1],na.rm = TRUE) / nrow(no_plasma_nc)

# NOW: run tmle 
whole = data.frame(whole)
whole$A = rep(0,77)
whole$id = 1:nrow(whole)

# rule_1_off: act>=128, plasma = 0
rule_1_off = whole[(whole[,2] >= 128 & whole[,3]==0),] 

# rule_1_on: act>=128, plasma > 0
rule_1_on = whole[(whole[,2] >= 128 & whole[,3]>0),] 

# rule_2_on: no act's or act<128, plasma = 0 
rule_2_on = whole[((whole[,2] < 128 | is.na(whole[,2])) & whole[,3]==0),] 

# rule_2_off: no act's or act<128, plasma > 0 
rule_2_off = whole[((whole[,2] < 128 | is.na(whole[,2])) & whole[,3]>0),] 

# off protocol = 0, on protocol = 1
A = rep(0,77)
A[rule_1_off$id] = 0
A[rule_1_on$id] = 1
A[rule_2_off$id] = 0
A[rule_2_on$id] = 1
whole$A = A


#unadjusted ate
t = whole[whole$A==1,]
# on/death
on_d = sum(t$final_out)/nrow(t)
# off/death
off_d = (nrow(t)-sum(t$final_out))/nrow(t)
# ate
on_d - off_d


# very raw tmle 
whole = whole[,-5]
whole = whole[-c(54,68,77),]

lib = c("SL.glm","SL.glmnet","SL.mean","SL.rpart","SL.xgboost","SL.ranger")
est_act = tmle(Y = whole$final_out,A = whole$A,W = whole[,4:9],Q.SL.library = lib,g.SL.library = lib)
est_act
