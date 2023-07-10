rm(list=ls())

#print(Sys.time())
#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)
library(BiocParallel)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(sandwich)
library(lmtest)
library(nnet)


cores=10

##Folder to save the outputs
output.folder="./results/"

# load data
dades=import("0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")

# regression functions
fastlm.fun <- function(x,y, z,x.var,y.var,sel) {

  tryCatch({

    y <- y
    x <-scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    prev=(table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))

    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(var.x=x.var,var.y=y.var,sel=sel,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data),message = NA)

  }, warning = function(w) {
    data.frame(var.x=x.var,var.y=y.var,sel=sel,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))

  }, error = function(e) {

    data.frame(var.x=x.var,var.y=y.var,sel=sel,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))

  })

}


#outcome
yi="BMI"
#covariates
covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")
sel="sbp_m_all"

#log1p transform
sbp_res=import("./results/res_fastlm.24SBP_multi_model.tsv")
sbp_res=sbp_res[which(sbp_res$q.value<0.05),]
noms=sbp_res$var.x

#run the models

system.time(res.lm1<-bplapply(noms,function(x){
  res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi,sel=sel)
  data.frame(res.lm1)
}, BPPARAM = MulticoreParam(cores)))
res.lm1 <- do.call(rbind, res.lm1)

sel="dbp_m_all"
dbp_res=import("./results/res_fastlm.24DBP_multi_model.tsv")
dbp_res=dbp_res[which(dbp_res$q.value<0.05),]
noms=dbp_res$var.x

#run the models

system.time(res.lm2<-bplapply(noms,function(x){
  res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi,sel=sel)
  data.frame(res.lm2)
}, BPPARAM = MulticoreParam(cores)))
res.lm2 <- do.call(rbind, res.lm2)

sel="syssd_all"
sysSD_res=import("./results/res_fastlm.sysSD_multi_model.tsv")
sysSD_res=sysSD_res[which(sysSD_res$q.value<0.05),]
noms=sysSD_res$var.x

#run the models

system.time(res.lm5<-bplapply(noms,function(x){
  res.lm5 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi,sel=sel)
  data.frame(res.lm5)
}, BPPARAM = MulticoreParam(cores)))
res.lm5 <- do.call(rbind, res.lm5)

#outcome

sel="diasd_all"
diaSD_res=import("./results/res_fastlm.diaSD_multi_model.tsv")
diaSD_res=diaSD_res[which(diaSD_res$q.value<0.05),]
noms=diaSD_res$var.x

#run the models

system.time(res.lm6<-bplapply(noms,function(x){
  res.lm6 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi,sel=sel)
  data.frame(res.lm6)
}, BPPARAM = MulticoreParam(cores)))
res.lm6 <- do.call(rbind, res.lm6)

#joint both phenotype analyses
res=data.frame(rbind(res.lm1,res.lm2,res.lm5,res.lm6),stringsAsFactors=F)
#estimate the fdr value
res$q.value <- p.adjust(res$p.value, method = "BH")

sbp=res[which(res$sel=="sbp_m_all"),]
sbp1 <- sbp[order(sbp$p.value),]

dbp=res[which(res$sel=="dbp_m_all"),]
dbp1 <- dbp[order(dbp$p.value),]

sbp_sd=res[which(res$sel=="syssd_all"),]
sbp_sd1 <- sbp_sd[order(sbp_sd$p.value),]

dbp_sd=res[which(res$sel=="diasd_all"),]
dbp_sd1 <- dbp_sd[order(dbp_sd$p.value),]


#save the result
write.table(sbp1,file=paste(output.folder,"/res_fastlm.BMI.24SBP_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(dbp1,file=paste(output.folder,"/res_fastlm.BMI.24DBP_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(sbp_sd1,file=paste(output.folder,"/res_fastlm.BMI.sysSD_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(dbp_sd1,file=paste(output.folder,"/res_fastlm.BMI.diaSD_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

