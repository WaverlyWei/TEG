###  Include simplest models that Houston Using ###

rm(list=ls())
setwd("~/Dropbox/hubbardlap/Mitch Cohen/Lucy/")
dataf=read.csv("acitLucy.csv")
#install.packages("BayesTree")
#library(BayesTree)
library(LogicReg)
library(ROCR)

### List of Outcomes
out=c("mt10rbcffp24h", "rbc_by_6h_bin", "ffp_by_6h_bin","mortalityat24h", "mortalityatdisch")	

stand.lab=c("hr0_inr","hr0_ptt")
standlabnmes=gsub("^.*?hr0_","",stand.lab)
#compared to the following groupings

#For TEG CK parameters:
#### Removed sample0h_ck_ly60 because almost no observations
teg.ck=c("sample0h_ck_sp","sample0h_ck_r","sample0h_ck_k","sample0h_ck_alpha","sample0h_ck_ma","sample0h_ck_g",
"sample0h_ck_ly30")
tegcknmes=gsub("^.*?sample0h_ck_","",teg.ck)

#For TEG CRT parameters:
teg.crt=c("sample0h_crt_act","sample0h_crt_sp","sample0h_crt_r","sample0h_crt_k","sample0h_crt_alpha","sample0h_crt_ma","sample0h_crt_g","sample0h_crt_ly30")
tegcrtnmes=gsub("^.*?sample0h_crt_","",teg.crt)

#For TEG FUNCTIONAL FIBRINOGEN parameters:
teg.ff=c("sample0h_ff_flev","sample0h_ff_ma")
tegffnmes=gsub("^.*?sample0h_ff_","",teg.ff)

#For ROTEM EXTEM parameters:
rotem.ex=c("sample0h_ex_ct","sample0h_ex_cft","sample0h_ex_alpha","sample0h_ex_a10","sample0h_ex_a20",
"sample0h_ex_mcf","sample0h_ex_ml")
rotemexnmes=gsub("^.*?sample0h_ex_","",rotem.ex)

#FOR ROTEM INTEM parameters:
rotem.intem=c("sample0h_in_ct","sample0h_in_cft","sample0h_in_alpha","sample0h_in_a10","sample0h_in_a20","sample0h_in_mcf","sample0h_in_ml")
rotemintemnmes=gsub("^.*?sample0h_in_","",rotem.intem)

#For ROTEM APTEM parameters:
rotem.aptem=c("sample0h_ap_ct","sample0h_ap_cft","sample0h_ap_alpha","sample0h_ap_a10","sample0h_ap_a20",
"sample0h_ap_mcf","sample0h_ap_ml")
rotemaptemnmes=gsub("^.*?sample0h_ap_","",rotem.aptem)

#For ROTEM FIBTEM parameters:
rotem.fibtem=c("sample0h_fib_ct","sample0h_fib_alpha","sample0h_fib_a10","sample0h_fib_a20","sample0h_fib_mcf",
"sample0h_fib_ml")
rotemfibtemnmes=gsub("^.*?sample0h_fib_","",rotem.fibtem)

rotem.all=c(rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem)
rotemnmes=c(rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes)
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
pred.sets.nmes=list(standlabnmes,tegcknmes,tegcrtnmes,tegffnmes,rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes,rotemnmes)

prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")
ns=length(pred.sets)
### Re-do outcomes
 Yt=dataf[,out[2]]
 n=length(Yt)
 Y=rep(NA,n)
 Y[Yt=="Yes"]=1
 Y[Yt=="No"]=0
 dataf[,out[2]]=Y

 Yt=dataf[,out[3]]
 n=length(Yt)
 Y=rep(NA,n)
 Y[Yt=="Yes"]=1
 Y[Yt=="No"]=0
 dataf[,out[3]]=Y

 Yt=dataf[,out[4]]
 n=length(Yt)
 Y=rep(NA,n)
 Y[Yt=="Dead"]=1
 Y[Yt=="Alive"]=0
 dataf[,out[4]]=Y

 Yt=dataf[,out[5]]
 n=length(Yt)
 Y=rep(NA,n)
 Y[Yt=="Dead"]=1
 Y[Yt=="Alive"]=0
 dataf[,out[5]]=Y

setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Lucy/Results/Logic Results/")

sum.na=function(x){sum(is.na(x))}

######### 
no=length(out)
for( i in 1:no) {
for(j in 1:ns) {
outc=out[i]
preds=pred.sets[[j]]
predn=prednmes[j]
predictorn=pred.sets.nmes[[j]]
cat(" i = ",i," out of  ",no,":  j = ",j," out of  ",ns,"\n")

Y=dataf[,outc]
include=is.na(Y)==F
dat.t=dataf[include,]

Y=dat.t[,outc]
X=dat.t[,preds]
Xna=apply(X,1,sum.na)
X=X[Xna==0,]
Y=Y[Xna==0]
filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

xp=dim(X)[2]
qt=apply(na.omit(X),2,quantile,probs=seq(0.1,0.9,0.2))
newX=NULL
coln=NULL
varn=predictorn
for(k in 1:xp) {
	xn=cut(X[,k],breaks=unique(c(min(X[,k])-0.1,qt[,k],max(X[,k]+0.1))))
	inds <- model.matrix(~ factor(xn) - 1)[,-1]
	nmes=colnames(inds)
	nc=nchar(nmes)
	nmes=paste(varn[k],substr(nmes,11,nc),sep="")
	coln=c(coln,nmes)
	newX=cbind(newX,inds)
}
colnames(newX)=coln
ncol=dim(newX)[2]
logicreg=logreg(Y, newX,select=1)
tit=paste("Predictors are ",predn," Outcome is ",outc,sep="")
pdf(paste(filepre,"pdf",sep="."))
plot.logregmodel2(logicreg$model,nms=coln,nny=0.25)
mtext(tit, side = 3, line = 1,cex=1.25)
dev.off()
}
}
