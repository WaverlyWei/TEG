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

## IMPUTATION
imp_W = amelia(W)
write.amelia(imp_W, file.stem = "imputed_data_set")
complete_W = read.csv("imputed_data_set1.csv")

# combine dataf with imputed W 
complete = data.frame(cbind(dataf,complete_W))
after_1_hr_death = complete[1:15,2272:2279]
# clean up outcome vars
grep("hours_to_death",nm)
## BUG: ?? imputed values disappear 
# ========= FIX THIS ===============
## MAYBE: dont's combine the above way?? 
## 
# remove early death & blank labs at 2&4 hr
after_1_hr_death = complete[(complete$hours_to_death>1),]
after_1_hr_death = complete[(complete$hours_to_death>1)&(!is.na(dataf$hr2_inr))&(!is.na(dataf$hr4_inr))&
                           (!is.na(dataf$hr2_ptt))&(!is.na(dataf$hr4_ptt)),]
# death =1  alive = 0 / only 0-6h
out = after_1_hr_death %>% dplyr::select(contains("mortalityat6h"))
## only 22 obs 
out$id = 1:nrow(out)
alive = out[out$mortalityat6h==0,]$id

## Targeting on MA
ma = after_1_hr_death %>% dplyr::select(contains("crt_ma"))
plt = after_1_hr_death%>% dplyr::select(contains("plt_units"))

# extract imputed_W
final_W = after_1_hr_death[,2272:2279]


# code treatment A
whole = cbind(final_out,final_ma,final_plt,final_W)
# comply
comply = whole[(whole[,2] < 55 & !is.na(whole[,2])),]
nrow(comply)

# comply + plt
get_plt = comply[(comply[,3]>0),]
nrow(get_plt)
# comply + no plt
no_plt = comply[(comply[,3]==0),]
nrow(no_plt)

# not comply
not_comply = whole[(whole[,2] >= 55 | is.na(whole[,2])),] 
nrow(not_comply)
# not comply + plt
get_plt_nc = not_comply[(not_comply[,3]>0),]
nrow(get_plt_nc)
# not comply + no plt
no_plt_nc = not_comply[(not_comply[,3]==0),]
nrow(no_plt_nc)


# death rate?? maybe not remove deaths?? CORRECT THIS 

sum(get_plt[,1],na.rm = TRUE) / nrow(get_plt)
sum(no_plt[,1],na.rm = TRUE) / nrow(no_plt)
sum(get_plt_nc[,1],na.rm = TRUE) / nrow(get_plt_nc)
sum(no_plt_nc[,1],na.rm = TRUE) / nrow(no_plt_nc)

# NOW: run tmle 
whole = data.frame(whole)
whole$A = rep(0,77)
whole$id = 1:nrow(whole)

# rule_1_off: MA<55, plt = 0
rule_1_off = whole[(whole[,2] < 55 & whole[,3]==0),] 

# rule_1_on: MA<55, plt>0
rule_1_on = whole[(whole[,2] < 55 & whole[,3]>0),] 

# rule_2_off:no MA's, plt > 0 
rule_2_off = whole[((whole[,2] >= 55 | is.na(whole[,2])) & whole[,3]>0),] 

# rule_2_on: no MA's, plt = 0 
rule_2_on = whole[((whole[,2] >= 55 | is.na(whole[,2])) & whole[,3]>0),] 

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
est_ma = tmle(Y = whole$final_out,A = whole$A,W = whole[,4:9],Q.SL.library = lib,g.SL.library = lib)
est_ma
