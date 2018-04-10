library(ggplot2)
library(wesanderson)
load("AUC_all.Rda")
load("auc_df.Rda")
df = auc_df_all
df$PredSet = rep("All",5)
###
group = factor(c(1,2,3,4,5))
df = cbind(df,group)
p <- ggplot(df,aes(df$PredSet,df$AUC,colour = ))
#p+ scale_y_discrete(breaks = seq(0,1,by = 0.001))
p + geom_linerange(aes(ymin = df$CI_low,ymax = df$CI_high))
p + geom_pointrange(aes(ymin = df$CI_low,ymax = df$CI_high))
#p + geom_crossbar(aes(ymin = df$CI_low,ymax = df$CI_high), width = 0.2)
p + geom_errorbar(aes(ymin = df$CI_low,ymax = df$CI_high, colour = group), width = 0.1) + theme(axis.text.y=element_blank())


p + geom_point() + geom_point(data = df_std,aes(df_std$PredSet,df_std$AUC,colour = group) )

# try 
try = rbind(df,df_std)
p <- ggplot(try,aes(try$PredSet,try$AUC,colour = group))
p + geom_errorbar(aes(ymin = try$CI_low,ymax = try$CI_high), width = 0.1)



##  run 
df_std = auc_df[auc_df$PredSet=="standLab",]
df_teg_ck = auc_df[auc_df$PredSet=="teg_ck",]
df_teg_crt = auc_df[auc_df$PredSet=="teg_crt",]
df_teg_ff = auc_df[auc_df$PredSet=="teg_ff",]
df_rotem_ex = auc_df[auc_df$PredSet=="rotem_ex",]
df_rotem_intem = auc_df[auc_df$PredSet=="rotem_intem",]
df_rotem_aptem = auc_df[auc_df$PredSet=="rotem_aptem",]
df_rotem_fibtem = auc_df[auc_df$PredSet=="rotem_fibtem",]
df_rotem_all = auc_df[auc_df$PredSet=="rotem_all",]

new_df = rbind(df,df_std,df_teg_ck,df_teg_crt,df_teg_ff,df_rotem_ex,df_rotem_intem,df_rotem_aptem,
               df_rotem_fibtem,df_rotem_all)
group = rep(c(1,2,3,4,5),10)
new_df = cbind(new_df,group)
p <- ggplot(new_df,aes(factor(new_df$PredSet),new_df$AUC)) + theme(axis.text.y=element_blank())
p + geom_errorbar(aes(ymin = new_df$CI_low,ymax = new_df$CI_high,colour = factor(group),width = 0.1))

