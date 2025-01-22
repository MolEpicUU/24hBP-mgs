print(Sys.time())

rm(list=ls())
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

cores=16

##Folder to save the outputs
output.folder="./results/"


# load data
dades=import("0_data/dades.pheno_malmo_2575_uppsala_195_validation.csv")
res1=import("results/res_fastlm.24SBP_multiBMI_model.tsv")
res2=import("results/res_fastlm.24DBP_multiBMI_model.tsv")

res1=res1[which(res1$p.value.infl<0.05 & sign(res1$estimate)==sign(res1$estimate.infl)),]
res2=res2[which(res2$p.value.infl<0.05 & sign(res2$estimate)==sign(res2$estimate.infl)),]

sel1=res1[which(res1$q.value<0.05),"var.x"]
sel2=res2[which(res2$q.value<0.05),"var.x"]

# regression functions
fastlm.fun <- function(x,y, z,x.var,y.var) {

  tryCatch({

    y <-y
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    prev= (table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data),message = NA)

  }, warning = function(w) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA,message=paste("Warning:", w$message))

  }, error = function(e) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))

  })

}
#variables to assess
noms=grep("hMGS.",names(dades),value=T)

#outcome
yi="bbps"
#covariates
covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI")

#run the models
print("##### time #####")
system.time(res.lm1<-bplapply(sel1,function(x){
  res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm1)
}, BPPARAM = MulticoreParam(cores)))
res.lm1 <- do.call(rbind, res.lm1)
res.lm1 =data.frame(res.lm1,stringsAsFactors=F)
#estimate the fdr value
res.lm1$q.value <- p.adjust(res.lm1$p.value, method = "BH")

#outcome
yi="bbpd"

#run the models
system.time(res.lm2<-bplapply(sel2,function(x){
  res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm2)
}, BPPARAM = MulticoreParam(cores)))
res.lm2 <- do.call(rbind, res.lm2)
res.lm2 =data.frame(res.lm2,stringsAsFactors=F)
#estimate the fdr value
res.lm2$q.value <- p.adjust(res.lm2$p.value, method = "BH")

sbp1 <- res.lm1[order(res.lm1$p.value),]
dbp1 <- res.lm2[order(res.lm2$p.value),]

#save the result
write.table(sbp1,file=paste(output.folder,"/res_fastlm.SBP_multiBMI_model_validation.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

write.table(dbp1,file=paste(output.folder,"/res_fastlm.DBP_multiBMI_model_validation.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

