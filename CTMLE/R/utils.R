getW <- function(){
  #W's
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
  #imp_W = amelia(W)
  #write.amelia(imp_W, file.stem = "imputed_data_set_alpha")
  complete_W <- read.csv("imputed_data_set_ma1.csv")
  return(cbind(complete_W,W_missing))
}


getA <- function(measure,product,threshold){
  ma <- dataf%>% dplyr::select(contains("crt_ma"))
  plt <- dataf %>% dplyr::select(contains("icu_0to6h_plt_units"))
  A <- rep(NA,nrow(dataf))
  #0-4h
  ma <- ma[,1:4]
  # form a temporary dataframe
  treatment <- cbind(ma,plt, A)
  
  comply <- row.names(treatment[apply(treatment[,1:4],1,
                                      function(x) any(x < 55,na.rm = TRUE)),])
  not_comply <- row.names(treatment[apply(treatment[,1:4],1,
                                          function(x) all(x >= 55,na.rm = TRUE)),])
  get_product <- row.names(treatment[(treatment[,5]>0 &!is.na(treatment[,5])),])
  no_product <- row.names(treatment[(treatment[,5]==0&!is.na(treatment[,5])),])
  
  
  
  # rule_1_on: ANY act>=128, plasma > 0
  rule.1.on <- as.numeric(unique(c(comply,get_product)))
  treatment[as.numeric(rule.1.on),ncol(treatment)] <- 1
  
  # rule_1_off: ANY act>=128, plasma = 0 
  rule.1.off <- as.numeric(unique(c(comply,no_product)))
  treatment[as.numeric(rule.1.off),ncol(treatment)] <- 0
  
  # rule_2_on: no act's or act<128, plasma = 0 
  rule.2.on <- as.numeric(unique(c(not_comply,no_product)))
  treatment[as.numeric(rule.2.on),ncol(treatment)] <- 1
  
  # rule_2_off: no act's or act<128, plasma > 0 
  rule.2.off <- as.numeric(unique(c(not_comply,get_product)))
  treatment[as.numeric(rule.2.off),ncol(treatment)] <- 0
  
  return(treatment$A)
}


# getY <- function(outcome,low,high){
#   bi.out <- ifelse(outcome>=low & outcome <= high,1,0)
#   return(bi.out)
# }

getY <- function(){
  Y <- rep(0,nrow(dataf))
  Y[which(dataf$icu_0to6h_blood_units>0 & dataf$icu_7to12h_blood_units==0)] <- 1
  return(Y)
}



getMort6 <- function(){
  # clean up outcome vars
  #m[grep("hours_to_death",nm)]
  # remove early death & blank labs at 2&4 hr
  dataf$id <- 1:nrow(dataf)
  # only remove early death FOR NOW
  after_1_hr_death = dataf[which(dataf$hours_to_death>1),]
  
  ## NOTE: if apply all the filter conditions, only 5 left 
  #after_1_hr_death = after_1_hr_death[((dataf$hours_to_death>1)&(!is.na(dataf$hr2_inr))&(!is.na(dataf$hr4_inr))&(!is.na(dataf$hr2_ptt))&(!is.na(dataf$hr4_ptt))),]
  
  index <- after_1_hr_death$id
  
  
  # death =1  alive = 0 / only 0-6h
  out = after_1_hr_death %>% dplyr::select(contains("mortalityat6h"))
  # 220 alive 57 dead
  #out$id = 1:nrow(out)
  return(list(out,index))
  #alive = out[out$mortalityat6h==0,]$id
}

getMort24 <- function(){
  # clean up outcome vars
  #m[grep("hours_to_death",nm)]
  # remove early death & blank labs at 2&4 hr
  dataf$id <- 1:nrow(dataf)
  # only remove early death FOR NOW
  after_1_hr_death = dataf[which(dataf$hours_to_death>1),]
  
  ## NOTE: if apply all the filter conditions, only 5 left 
  #after_1_hr_death = after_1_hr_death[((dataf$hours_to_death>1)&(!is.na(dataf$hr2_inr))&(!is.na(dataf$hr4_inr))&(!is.na(dataf$hr2_ptt))&(!is.na(dataf$hr4_ptt))),]
  
  index <- after_1_hr_death$id
  
  
  # death =1  alive = 0 / only 0-6h
  out = after_1_hr_death %>% dplyr::select(contains("mortalityat24h"))
  # 220 alive 57 dead
  #out$id = 1:nrow(out)
  return(list(out,index))
  #alive = out[out$mortalityat6h==0,]$id
}


getMortDisch <- function(){
  # clean up outcome vars
  #m[grep("hours_to_death",nm)]
  # remove early death & blank labs at 2&4 hr
  dataf$id <- 1:nrow(dataf)
  # only remove early death FOR NOW
  after_1_hr_death = dataf[which(dataf$hours_to_death>1),]
  
  ## NOTE: if apply all the filter conditions, only 5 left 
  #after_1_hr_death = after_1_hr_death[((dataf$hours_to_death>1)&(!is.na(dataf$hr2_inr))&(!is.na(dataf$hr4_inr))&(!is.na(dataf$hr2_ptt))&(!is.na(dataf$hr4_ptt))),]
  
  index <- after_1_hr_death$id
  
  
  # death =1  alive = 0 / only 0-6h
  out = after_1_hr_death %>% dplyr::select(contains("mortalityatdisch"))
  # 220 alive 57 dead
  #out$id = 1:nrow(out)
  return(list(out,index))
  #alive = out[out$mortalityat6h==0,]$id
}





UATE <- function(df){
  
}