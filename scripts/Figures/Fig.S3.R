
rm(list=ls())

library(rio)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggforestplot)
library(ggeasy)
#output folder

sbp.reg.mod1=sbp.reg.mod1[which(sbp.reg.mod1$q.value<0.05),]
dbp.reg.mod1=dbp.reg.mod1[which(dbp.reg.mod1$q.value<0.05),]
sbp.sd.reg.mod1=sbp.sd.reg.mod1[which(sbp.sd.reg.mod1$q.value<0.05),]
dbp.sd.reg.mod1=dbp.sd.reg.mod1[which(dbp.sd.reg.mod1$q.value<0.05),]

sbp.reg.mod1$var.x<- gsub("____", " (", sbp.reg.mod1$var.x)
sbp.reg.mod1$var.x<- gsub("_", " ", sbp.reg.mod1$var.x) 
sbp.reg.mod1$var.x<- gsub("$", ")", sbp.reg.mod1$var.x)


dbp.reg.mod1$var.x<- gsub("____", " (", dbp.reg.mod1$var.x)
dbp.reg.mod1$var.x<- gsub("_", " ", dbp.reg.mod1$var.x) 
dbp.reg.mod1$var.x<- gsub("$", ")", dbp.reg.mod1$var.x)

sbp.sd.reg.mod1$var.x<- gsub("____", " (", sbp.sd.reg.mod1$var.x)
sbp.sd.reg.mod1$var.x<- gsub("_", " ", sbp.sd.reg.mod1$var.x) 
sbp.sd.reg.mod1$var.x<- gsub("$", ")", sbp.sd.reg.mod1$var.x)

dbp.sd.reg.mod1$var.x<- gsub("____", " (", dbp.sd.reg.mod1$var.x)
dbp.sd.reg.mod1$var.x<- gsub("_", " ", dbp.sd.reg.mod1$var.x) 
dbp.sd.reg.mod1$var.x<- gsub("$", ")", dbp.sd.reg.mod1$var.x)




#sbp
sbp.reg.mod1$model="main_regression_model"
sbp.reg.s1$model="sensitivity_model"
sbp_all=merge(sbp.reg.mod1,sbp.reg.s1, by="var.x")

#dbp
dbp.reg.mod1$model="main_regression_model"
dbp.reg.s1$model="sensitivity_model"
dbp_all=merge(dbp.reg.mod1,dbp.reg.s1, by="var.x")


#sbp-sd
sbp.sd.reg.mod1$model="main_regression_model"
sbp.sd.reg.s1$model="sensitivity_model"
sbp_sd_all=merge(sbp.sd.reg.mod1,sbp.sd.reg.s1, by="var.x")

#dbp-sd
dbp.sd.reg.mod1$model="main_regression_model"
dbp.sd.reg.s1$model="sensitivity_model"
dbp_sd_all=merge(dbp.sd.reg.mod1,dbp.sd.reg.s1, by="var.x")

#################
s1=ggplot(sbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=3886)")+ylab("Estimates in sensitivity 1")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  ggtitle("SBP(sensitivity 1)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-50), label.y = 20)
#aes(label = ..r.label..),
s1

s2=ggplot(dbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=3886)")+ylab("Estimates in sensitivity 1")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+

  ggtitle("DBP(sensitivity 1)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-35), label.y = 4.9)
#(aes(label = ..r.label..),
s2

s3=ggplot(sbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=3886)")+ylab("Estimates in sensitivity 1")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
 
ggtitle("SBP-SD(sensitivity 1)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-8), label.y = 15)
#aes(label = ..r.label..),
s3

s4=ggplot(dbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=3886)")+ylab("Estimates in sensitivity 1")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
ggtitle("DBP-SD(sensitivity 1)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-10), label.y = 4)
s4


p1=ggarrange(s1,s2,s3,s4,ncol=2, nrow=2,common.legend = TRUE, legend="top")
p1


##################################################
##################################################
##################################################


#import results

sbp.reg.mod1=sbp.reg.mod1[which(sbp.reg.mod1$q.value<0.05),]
dbp.reg.mod1=dbp.reg.mod1[which(dbp.reg.mod1$q.value<0.05),]
sbp.sd.reg.mod1=sbp.sd.reg.mod1[which(sbp.sd.reg.mod1$q.value<0.05),]
dbp.sd.reg.mod1=dbp.sd.reg.mod1[which(dbp.sd.reg.mod1$q.value<0.05),]

sbp.reg.mod1$q.value=-log10(sbp.reg.mod1$q.value)
dbp.reg.mod1$q.value=-log10(dbp.reg.mod1$q.value)
sbp.sd.reg.mod1$q.value=-log10(sbp.sd.reg.mod1$q.value)
dbp.sd.reg.mod1$q.value=-log10(dbp.sd.reg.mod1$q.value)

sbp.reg.s1$q.value=-log10(sbp.reg.s1$q.value)
dbp.reg.s1$q.value=-log10(dbp.reg.s1$q.value)
sbp.sd.reg.s1$q.value=-log10(sbp.sd.reg.s1$q.value)
dbp.sd.reg.s1$q.value=-log10(dbp.sd.reg.s1$q.value)

#adjust name
sbp.reg.mod1$var.x<- gsub("____", " (", sbp.reg.mod1$var.x)
sbp.reg.mod1$var.x<- gsub("_", " ", sbp.reg.mod1$var.x) 
sbp.reg.mod1$var.x<- gsub("$", ")", sbp.reg.mod1$var.x)


dbp.reg.mod1$var.x<- gsub("____", " (", dbp.reg.mod1$var.x)
dbp.reg.mod1$var.x<- gsub("_", " ", dbp.reg.mod1$var.x) 
dbp.reg.mod1$var.x<- gsub("$", ")", dbp.reg.mod1$var.x)

sbp.sd.reg.mod1$var.x<- gsub("____", " (", sbp.sd.reg.mod1$var.x)
sbp.sd.reg.mod1$var.x<- gsub("_", " ", sbp.sd.reg.mod1$var.x) 
sbp.sd.reg.mod1$var.x<- gsub("$", ")", sbp.sd.reg.mod1$var.x)

dbp.sd.reg.mod1$var.x<- gsub("____", " (", dbp.sd.reg.mod1$var.x)
dbp.sd.reg.mod1$var.x<- gsub("_", " ", dbp.sd.reg.mod1$var.x) 
dbp.sd.reg.mod1$var.x<- gsub("$", ")", dbp.sd.reg.mod1$var.x)


#sbp
sbp.reg.mod1$model="main_regression_model"
sbp.reg.s1$model="sensitivity_model"
sbp_all=merge(sbp.reg.mod1,sbp.reg.s1, by="var.x")

#dbp
dbp.reg.mod1$model="main_regression_model"
dbp.reg.s1$model="sensitivity_model"
dbp_all=merge(dbp.reg.mod1,dbp.reg.s1, by="var.x")


#sbp-sd
sbp.sd.reg.mod1$model="main_regression_model"
sbp.sd.reg.s1$model="sensitivity_model"
sbp_sd_all=merge(sbp.sd.reg.mod1,sbp.sd.reg.s1, by="var.x")

#dbp-sd
dbp.sd.reg.mod1$model="main_regression_model"
dbp.sd.reg.s1$model="sensitivity_model"
dbp_sd_all=merge(dbp.sd.reg.mod1,dbp.sd.reg.s1, by="var.x")

#################
s1=ggplot(sbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4253)")+ylab("Estimates in sensitivity 2")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
 ggtitle("SBP(sensitivity 2)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-50), label.y = 20)
#aes(label = ..r.label..),
s1

s2=ggplot(dbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4253)")+ylab("Estimates in sensitivity 2")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("DBP(sensitivity 2)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-35), label.y = 5)
#(aes(label = ..r.label..),
s2

s3=ggplot(sbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4253)")+ylab("Estimates in sensitivity 2")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("SBP-SD(sensitivity 2)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-8), label.y = 15)
#aes(label = ..r.label..),
s3

s4=ggplot(dbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4253)")+ylab("Estimates in sensitivity 2")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("DBP-SD(sensitivity 2)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-10), label.y = 4)
#(aes(label = ..r.label..),
s4


p2=ggarrange(s1,s2,s3,s4,ncol=2, nrow=2,common.legend = TRUE, legend="top")
p2

#####################################################
####################################################
####################################################
####################################################
#import results
sbp.reg.mod1=sbp.reg.mod1[which(sbp.reg.mod1$q.value<0.05),]
dbp.reg.mod1=dbp.reg.mod1[which(dbp.reg.mod1$q.value<0.05),]
sbp.sd.reg.mod1=sbp.sd.reg.mod1[which(sbp.sd.reg.mod1$q.value<0.05),]
dbp.sd.reg.mod1=dbp.sd.reg.mod1[which(dbp.sd.reg.mod1$q.value<0.05),]
ㄆ
sbp.reg.mod1$q.value=-log10(sbp.reg.mod1$q.value)
dbp.reg.mod1$q.value=-log10(dbp.reg.mod1$q.value)
sbp.sd.reg.mod1$q.value=-log10(sbp.sd.reg.mod1$q.value)
dbp.sd.reg.mod1$q.value=-log10(dbp.sd.reg.mod1$q.value)


sbp.reg.s1$q.value=-log10(sbp.reg.s1$q.value)
dbp.reg.s1$q.value=-log10(dbp.reg.s1$q.value)
sbp.sd.reg.s1$q.value=-log10(sbp.sd.reg.s1$q.value)
dbp.sd.reg.s1$q.value=-log10(dbp.sd.reg.s1$q.value)

#adjust name
sbp.reg.mod1$var.x<- gsub("____", " (", sbp.reg.mod1$var.x)
sbp.reg.mod1$var.x<- gsub("_", " ", sbp.reg.mod1$var.x) 
sbp.reg.mod1$var.x<- gsub("$", ")", sbp.reg.mod1$var.x)


dbp.reg.mod1$var.x<- gsub("____", " (", dbp.reg.mod1$var.x)
dbp.reg.mod1$var.x<- gsub("_", " ", dbp.reg.mod1$var.x) 
dbp.reg.mod1$var.x<- gsub("$", ")", dbp.reg.mod1$var.x)

sbp.sd.reg.mod1$var.x<- gsub("____", " (", sbp.sd.reg.mod1$var.x)
sbp.sd.reg.mod1$var.x<- gsub("_", " ", sbp.sd.reg.mod1$var.x) 
sbp.sd.reg.mod1$var.x<- gsub("$", ")", sbp.sd.reg.mod1$var.x)

dbp.sd.reg.mod1$var.x<- gsub("____", " (", dbp.sd.reg.mod1$var.x)
dbp.sd.reg.mod1$var.x<- gsub("_", " ", dbp.sd.reg.mod1$var.x) 
dbp.sd.reg.mod1$var.x<- gsub("$", ")", dbp.sd.reg.mod1$var.x)



#sbp
sbp.reg.mod1$model="main_regression_model"
sbp.reg.s1$model="sensitivity_model"
sbp_all=merge(sbp.reg.mod1,sbp.reg.s1, by="var.x")

#dbp
dbp.reg.mod1$model="main_regression_model"
dbp.reg.s1$model="sensitivity_model"
dbp_all=merge(dbp.reg.mod1,dbp.reg.s1, by="var.x")


#sbp-sd
sbp.sd.reg.mod1$model="main_regression_model"
sbp.sd.reg.s1$model="sensitivity_model"
sbp_sd_all=merge(sbp.sd.reg.mod1,sbp.sd.reg.s1, by="var.x")

#dbp-sd
dbp.sd.reg.mod1$model="main_regression_model"
dbp.sd.reg.s1$model="sensitivity_model"
dbp_sd_all=merge(dbp.sd.reg.mod1,dbp.sd.reg.s1, by="var.x")

#################
s1=ggplot(sbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4200)")+ylab("Estimates in sensitivity 3")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("SBP(sensitivity 3)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-50), label.y = 20)
#aes(label = ..r.label..),
s1

s2=ggplot(dbp_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4200)")+ylab("Estimates in sensitivity 3")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("DBP(sensitivity 3)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-35), label.y = 5)
#(aes(label = ..r.label..),
s2

s3=ggplot(sbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4200)")+ylab("Estimates in sensitivity 3")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("SBP-SD(sensitivity 3)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor( digits=2,method = "spearman", cor.coef.name="rho", label.x =(-8), label.y = 15)
#aes(label = ..r.label..),
s3

s4=ggplot(dbp_sd_all,aes(x=estimate.x,y=estimate.y,label=var.x))+
  geom_point(color="cornflowerblue")+
  xlab("Estimates in main model(n=4200)")+ylab("Estimates in sensitivity 3")+
  theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )+
  #scale_y_continuous(limits = c(0, 5))+
  #scale_x_continuous(limits = c(0, 6))+
  
  # xlab(paste("estimate (all; n=", substring(aa$n.x,1,1),",",substring(aa$n.x,2,nchar(aa$n.x)),")",sep=""))+
  # ylab(paste("estimate (no crohn disease; n=", substring(aa$n.y,1,1),",",substring(aa$n.y,2,nchar(aa$n.y)),")",sep=""))+
  # geom_label_repel(data=head(aa[order(aa$diff.est,decreasing = T),],dif.n),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  force=70)+
  # 
ggtitle("DBP-SD(sensitivity 3)") +
  ggeasy::easy_center_title()+
  geom_smooth(method="lm", se=F, colour = "gray" )+stat_cor(digits=2,method = "spearman", cor.coef.name="rho", label.x =(-10), label.y = 4)
#(aes(label = ..r.label..),
s4


p3=ggarrange(s1,s2,s3,s4,ncol=2, nrow=2,common.legend = TRUE, legend="top")
p3
####################
p5=ggarrange(p1,p2,p3, ncol=1, nrow=3,common.legend = TRUE, legend="top")
p5



