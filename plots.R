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
p = p + geom_errorbar(aes(ymin = df$CI_low,ymax = df$CI_high, colour = group), width = 0.1) + theme(axis.text.y=element_blank())
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
outcome = rep(c(1,2,3,4,5),10)
#outcome = factor(outcome)

#levels(new_df$outcome)[levels(new_df$outcome)=="1"] <- "Control"
#levels(pg$group)[levels(pg$group)=="trt1"] <- "Treatment 1"
#levels(pg$group)[levels(pg$group)=="trt2"] <- "Treatment 2"
#names(pg)[names(pg)=="group"]  <- "Experimental Condition"



new_df = cbind(new_df,outcome)
#levels(new_df$outcome)[levels(new_df$outcome)==1] <- "Control"

new_df$AUC = as.numeric(new_df$AUC)
new_df$CI_low = as.numeric(new_df$CI_low)
new_df$CI_high = as.numeric(new_df$CI_high)
new_df = data.frame(new_df, stringsAsFactors = FALSE)

p = ggplot(new_df,aes(factor(new_df$PredSet),new_df$AUC))  
# theme(axis.text.y=element_blank()) # hide y axis
p = p + geom_errorbar(aes(ymin = new_df$CI_low,ymax = new_df$CI_high,colour = factor(outcome),width = 0.15))+
  ylab("AUC")+xlab("PredSet")+ggtitle("AUC of All PredSets")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
p = p + geom_hline(yintercept = 0.5, colour = "red") + ylim(0,1) + theme(axis.text.x = element_text(face="bold", size=8, angle=45))
p = p + theme(plot.title = element_text(hjust = 0.5)) 






# heat map
h = ggplot(new_df, aes(PredSet, Outcome)) +
  geom_tile(aes(fill = AUC), color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Outcome ") +
  xlab("PredSet") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "AUC")
