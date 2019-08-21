rm(list=ls())
library(MASS)
library(SuperLearner)
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
library(speedglm)
library(sl3)
library(origami)
library(devtools)
library(optxx)

set.seed(12313197)
options(digits = 4)
#setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Trauma and Coagulation (White Space Conflict)/ACIT/Data")
setwd('/Users/waverlywei/Desktop/Lucy project')
dataf=read.dta13("ACIT_6Jan17_deidentified1.dta")
plt0to6=as.numeric(dataf$icu_0to6h_plt_units>0)
bloodUnit0to6=as.numeric(dataf$icu_0to6h_blood_units>0)
ffp0to6=as.numeric(dataf$icu_0to6h_ffp_units>0)
dataf=data.frame(dataf,plt0to6,bloodUnit0to6,ffp0to6)
rm(plt0to6,bloodUnit0to6,ffp0to6)

## ======CODE W's ====== ##

# locate covariates 
nm = names(dataf)
# check position 
nm[grep("iss",nm)]
nm[grep("hr0_basedefexc",nm)]
grep("male",nm)
grep("age",nm)
# mechanism vars
grep("blunt",nm)

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
# ADD Mechtype / blunt
mech = as.character(dataf[,514])
mech[mech=="Blunt"] = 0
mech[mech=="Penetrating"] = 1
mech = as.numeric(mech)
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
W = cbind(iss,basexc,gender,age,mech,sbp,Hr,race)
# create missing indicators
W_missing = W
# missing = 1
W_missing = apply(W_missing,2,is.na)
W_missing = apply(W_missing,2,as.numeric)
colnames(W_missing) = sapply(colnames(W_missing),function(x) paste0("missing_",x))

## IMPUTATION
imp_W = amelia(W)
write.amelia(imp_W, file.stem = "imputed_data_set_ma")
complete_W = read.csv("imputed_data_set_ma1.csv")

## final W
final_W = cbind(complete_W,W_missing)



## ======CODE V's ====== #
# ffp as product 

V = dataf %>% select(hr0_inr, hr0_ptt, sample0h_ex_ct, sample0h_crt_r, 
                      sample0h_ck_r, sample0h_crt_act, sample0h_ex_cft, 
                      sample0h_crt_k, sample0h_ck_k)

# create missing indicators
V_missing = V
# missing = 1
V_missing <- apply(V_missing,2,is.na)
V_missing <- apply(V_missing,2,as.numeric)
colnames(V_missing) <- sapply(colnames(V_missing),function(x) paste0("missing_",x))

## IMPUTATION
imp_V <- amelia(V)
write.amelia(imp_V, file.stem ="imputed_data_set_V")
complete_V <- read.csv("imputed_data_set_V1.csv")

## final W
final_V <- cbind(complete_V,V_missing)

## code Vnorm 
# normal = 1   
Vnorm <- final_V
Vnorm  <- Vnorm %>% 
  mutate(sample0h_ex_ct=as.numeric((sample0h_ex_ct>=38) & (sample0h_ex_ct<=79))) %>% 
  mutate(sample0h_crt_r=as.numeric((sample0h_crt_r>=0) & (sample0h_crt_r<=1))) %>% 
  mutate(sample0h_ck_r=as.numeric((sample0h_ck_r>=2) & (sample0h_ck_r<=8))) %>% 
  mutate(sample0h_crt_act=as.numeric((sample0h_crt_act>=76) & (sample0h_crt_act<=110))) %>%
  mutate(sample0h_ex_cft=as.numeric((sample0h_ex_cft>=34) & (sample0h_ex_cft<=159))) %>% 
  mutate(sample0h_crt_k=as.numeric((sample0h_crt_k>=1) & (sample0h_crt_k<=2))) %>% 
  mutate(sample0h_ck_k=as.numeric((sample0h_ck_k>=1) & (sample0h_ck_k<=3))) 


## CODE A
A <- dataf %>% 
  select(icu_0to6h_cryo_units)%>% 
  mutate(icu_0to6h_cryo_units = as.numeric(icu_0to6h_cryo_units>0))

## CODE Y
# clean up outcome vars
grep("hours_to_death",nm)
# remove early death & blank labs at 2&4 hr
dataf$id = 1:nrow(dataf)
# only remove early death FOR NOW
after_1_hr_death = dataf[which(dataf$hours_to_death>1),]

## NOTE: if apply all the filter conditions, only 5 left 
#after_1_hr_death = after_1_hr_death[((dataf$hours_to_death>1)&(!is.na(dataf$hr2_inr))&(!is.na(dataf$hr4_inr))&(!is.na(dataf$hr2_ptt))&(!is.na(dataf$hr4_ptt))),]

index = after_1_hr_death$id

# death =1  alive = 0 / only 0-6h
out = after_1_hr_death %>% dplyr::select(contains("mortalityat6h"))
out$id = 1:nrow(out)
# 220 alive 57 dead
alive = out[out$mortalityat6h==0,]$id



all <- cbind(final_W[index,], Vnorm[index,],A[index,],out)
W_nam <- c("W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12","W13","W14","W15","W16")
V_nam <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18")
all <- all[,c(-1,-18,-39)]
names(all) <- c(W_nam,V_nam,"A","Y")
all <- all %>% drop_na()



opt <- tmleopttx(all,Vnodes = grep("^V",names(all),value = T))

rule_ffp <- opt$rule
tmle_ffp <- data.frame(cbind(opt$tmlePsi,opt$tmleSD,opt$tmleCI))

xtable(head(rule_ffp))
xtable(tmle_ffp)

# CHECK FOR NA's
apply(all,2,function(x) any(is.na(x)))

