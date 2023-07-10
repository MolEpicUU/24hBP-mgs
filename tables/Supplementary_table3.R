print(Sys.time())

rm(list=ls())

library(rio)

info <- import("raw/HG3.A.7_tax.xlsx")
colnames(info)[1:3] <- c("id", "name", "level")

sbp=import("results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp=import("results/res_fastlm.24DBP_multiBMI_model.tsv")
sbp_sd=import("results/res_fastlm.sysSD_multiBMI_model.tsv")
dbp_sd=import("results/res_fastlm.diaSD_multiBMI_model.tsv")

sel=unique(c(sbp$var.x,dbp$var.x,sbp_sd$var.x,dbp_sd$var.x))
sel=substring(sel,nchar(sel)-8, nchar(sel))

sbp$prev=format(round(sbp$prev,2),nsmall=2)
sbp$id=substring(sbp$var.x,nchar(sbp$var.x)-8, nchar(sbp$var.x))
sbp$estimate=format(round(sbp[,"estimate"],2),nsmall=2)
sbp$lower=format(round(sbp[,"lower"],2),nsmall=2)
sbp$upper=format(round(sbp[,"upper"],2),nsmall=2)
sbp$estimate.infl=format(round(sbp[,"estimate.infl"],2),nsmall=2)
sbp$lower.infl=format(round(sbp[,"lower.infl"],2),nsmall=2)
sbp$upper.infl=format(round(sbp[,"upper.infl"],2),nsmall=2)


dbp$id=substring(dbp$var.x,nchar(dbp$var.x)-8, nchar(dbp$var.x))
dbp$estimate=format(round(dbp[,"estimate"],2),nsmall=2)
dbp$lower=format(round(dbp[,"lower"],2),nsmall=2)
dbp$upper=format(round(dbp[,"upper"],2),nsmall=2)
dbp$estimate.infl=format(round(dbp[,"estimate.infl"],2),nsmall=2)
dbp$lower.infl=format(round(dbp[,"lower.infl"],2),nsmall=2)
dbp$upper.infl=format(round(dbp[,"upper.infl"],2),nsmall=2)

sbp_sd$id=substring(sbp_sd$var.x,nchar(sbp_sd$var.x)-8, nchar(sbp_sd$var.x))
sbp_sd$estimate=format(round(sbp_sd[,"estimate"],2),nsmall=2)
sbp_sd$lower=format(round(sbp_sd[,"lower"],2),nsmall=2)
sbp_sd$upper=format(round(sbp_sd[,"upper"],2),nsmall=2)
sbp_sd$estimate.infl=format(round(sbp_sd[,"estimate.infl"],2),nsmall=2)
sbp_sd$lower.infl=format(round(sbp_sd[,"lower.infl"],2),nsmall=2)
sbp_sd$upper.infl=format(round(sbp_sd[,"upper.infl"],2),nsmall=2)

dbp_sd$id=substring(dbp_sd$var.x,nchar(dbp_sd$var.x)-8, nchar(dbp_sd$var.x))
dbp_sd$estimate=format(round(dbp_sd[,"estimate"],2),nsmall=2)
dbp_sd$lower=format(round(dbp_sd[,"lower"],2),nsmall=2)
dbp_sd$upper=format(round(dbp_sd[,"upper"],2),nsmall=2)
dbp_sd$estimate.infl=format(round(dbp_sd[,"estimate.infl"],2),nsmall=2)
dbp_sd$lower.infl=format(round(dbp_sd[,"lower.infl"],2),nsmall=2)
dbp_sd$upper.infl=format(round(dbp_sd[,"upper.infl"],2),nsmall=2)

info=info[which(info$id%in%sel==T),]

basic=merge(info[,c("name","id")],sbp[,c("id","prev","n","estimate","lower","upper","p.value","q.value","n.infl","estimate.infl","lower.infl","upper.infl","p.value.infl","q.value.infl")],by="id",all.x=T)
basic=merge(basic,dbp[,c("id","n","estimate","lower","upper","p.value","q.value","n.infl","estimate.infl","lower.infl","upper.infl","p.value.infl","q.value.infl")],by="id",all.x=T)
basic=merge(basic,sbp_sd[,c("id","n","estimate","lower","upper","p.value","q.value","n.infl","estimate.infl","lower.infl","upper.infl","p.value.infl","q.value.infl")],by="id",all.x=T)
basic=merge(basic,dbp_sd[,c("id","n","estimate","lower","upper","p.value","q.value","n.infl","estimate.infl","lower.infl","upper.infl","p.value.infl","q.value.infl")],by="id",all.x=T)

basic=basic[order(basic[,8]),]

names(basic)[names(basic)%in%grep("value",names(basic),value=T)]=paste(grep("value",names(basic),value=T),1:length(names(basic)[names(basic)%in%grep("value",names(basic),value=T)]),sep="")

for (i in grep("value",names(basic),value=T))
{
  basic[,i]=ifelse(basic[,i]<0.001,formatC(basic[,i],format="e",digits=2),round(basic[,i],digit=3))
}

names(basic)=c("id","name","prevalence","n.sbp","estimate.sbp","lower.sbp","upper.sbp","p.value.sbp","q.value.sbp",
"n.infl.sbp","estimate.infl.sbp","lower.infl.sbp","upper.infl.sbp","p.value.infl.sbp","q.value.infl.sbp",
"n.dbp","estimate.dbp","lower.dbp","upper.dbp","p.value.dbp","q.value.dbp",
"n.infl.dbp","estimate.infl.dbp","lower.infl.dbp","upper.infl.dbp","p.value.infl.dbp","q.value.infl.dbp",
"n.sbp_sd","estimate.sbp_sd","lower.sbp_sd","upper.sbp_sd","p.value.dbp_sd","q.value.dbp",
"n.infl.sbp_sd","estimate.infl.sbp_sd","lower.infl.sbp_sd","upper.infl.sbp_sd","p.value.infl.sbp_sd","q.value.infl.sbp_sd",
"n.dbp_sd","estimate.dbp_sd","lower.dbp_sd","upper.dbp_sd","p.value.dbp_sd","q.value.dbp",
"n.infl.dbp_sd","estimate.infl.dbp_sd","lower.infl.dbp_sd","upper.infl.dbp_sd","p.value.infl.dbp_sd","q.value.infl.dbp_sd")

write.table(basic,file="results/tables/Supplementary_table3.tsv",col.names=T,row.names=F,sep="\t")
print(Sys.time())


