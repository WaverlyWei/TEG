## libraries 
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
library(rpgm)

## data
set.seed(12313197)
options(digits = 4)
#setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Trauma and Coagulation (White Space Conflict)/ACIT/Data")
setwd('/Users/weilinqing/Desktop/Lucy project')
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

# clean up outcome vars
grep("hours_to_death",nm)
# remove early death & blank labs at 2&4 hr
dataf$id = 1:nrow(dataf)
# only remove early death FOR NOW
after_1_hr_death = dataf[which(dataf$hours_to_death>1),]

index = after_1_hr_death$id

# death =1  alive = 0 / only 0-6h
out = after_1_hr_death %>% dplyr::select(contains("mortalityat6h"))
out$id = 1:nrow(out)
# 220 alive 57 dead
alive = out[out$mortalityat6h==0,]$id

# clean out 
out = out[,1]

# extract imputed_W
final_W = cbind(complete_W,W_missing)[index,]

# create library
lib = c("SL.glm","SL.glmnet","SL.mean","SL.rpart","SL.xgboost","SL.ranger")

## write estInit function 
# Use this for Initial estimates of Q and g 
estInit <- function(Y, X, folds = 10,library){
  all<-cbind.data.frame(Y,X)
  names(all) = c("out",names(X))
  print(names(all))
  
  #Create sl3 task:
  task <- sl3::make_sl3_Task(all, covariates = names(X), outcome = "out", folds=folds)
  glm_learner <- Lrnr_glm$new()
  glmnet_learner = Lrnr_glmnet$new()
  rpart_learner = Lrnr_rpart$new()
  mean_learner = Lrnr_mean$new()
  learner_stack = Stack$new(glm_learner,glmnet_learner,rpart_learner,mean_learner)
  stack_fit <- learner_stack$train(task)
   preds = stack_fit$predict()
  return(preds)
  #out <- list(valY = sl3_fit$pred, SL.library = library, folds = folds,
             # fullFit = sl3_fit$sl.fit, cvFit=sl3_fit$cv.fit)
  
}

# test estInit 
init_res = estInit(out,final_W,lib)


## write ruletmle 
ruletmle <-function(obsA, obsY, pA1, Q0W, Q1W, ruleA, maxIter=1000,Qbounds=c(1e-4,1-1e-4)){
  
}

## ================TRY PACKAGE =============== ##
library(optxx)
