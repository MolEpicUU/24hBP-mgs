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
fastlm.fun <- function(x,y, z,x.var,y.var) {
  
  tryCatch({
    
    y <- y
    x <-scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    prev=(table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))
    
    df<-dfbetas(fit)
    nb=which.is.max(abs(df[,2]))
    fit2 <- lm(y ~ ., data[-nb,])
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    coef2 <- summary(fit2)$coefficients
    ci2 <- confint(fit2)
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data), estimate.infl = coef2["x", 1], lower.infl = ci2["x", 1], upper.infl = ci2["x", 2], se.infl = coef2["x", 2], p.value.infl = coef2["x", 4], n.infl = nrow(data[-nb,]),message = NA)
    
  }, warning = function(w) {
    
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA,estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA,estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA, message = paste("Error:", e$message))
    
  })
  
}


#outcome
yi="bbps"
#covariates
covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")


#log1p transform
sbp_res=import("./results/res_fastlm.24SBP_crude_model.tsv")
sbp_res=sbp_res[which(sbp_res$q.value<0.05),]
noms=sbp_res$var.x

#run the models

system.time(res.lm1<-bplapply(noms,function(x){
  res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm1)
}, BPPARAM = MulticoreParam(cores)))
res.lm1 <- do.call(rbind, res.lm1)

#outcome
yi="bbpd"

dbp_res=import("./results/res_fastlm.24DBP_crude_model.tsv")
dbp_res=dbp_res[which(dbp_res$q.value<0.05),]
noms=dbp_res$var.x

#run the models

system.time(res.lm2<-bplapply(noms,function(x){
  res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
  data.frame(res.lm2)
}, BPPARAM = MulticoreParam(cores)))
res.lm2 <- do.call(rbind, res.lm2)

#joint both phenotype analyses
res=data.frame(rbind(res.lm1,res.lm2),stringsAsFactors=F)
#estimate the fdr value
res$q.value <- p.adjust(res$p.value, method = "BH")
res$q.value.infl<-p.adjust(res$p.value.infl,method="BH")

sbp=res[which(res$var.y=="bbps"),]
sbp1 <- sbp[order(sbp$p.value),]

dbp=res[which(res$var.y=="bbpd"),]
dbp1 <- dbp[order(dbp$p.value),]


#save the result
write.table(sbp1,file=paste(output.folder,"/res_fastlm.office_SBP_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")
write.table(dbp1,file=paste(output.folder,"/res_fastlm.office_DBP_multi_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")


