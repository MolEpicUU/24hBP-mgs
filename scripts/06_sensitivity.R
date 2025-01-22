rm(list=ls())
# load libraries and set options

library(rio)
library(BiocParallel)

cores <- 16
set.seed(1)

# import and clean data

data<- import("0_data/dades.pheno_malmo_758_uppsala_2937_new.csv")
info <- import("raw/scapis_metagenomics_mgs_annotations_v1.0.tsv")
sbp1 <- import("results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp1 <- import("results/res_fastlm.24DBP_multiBMI_model.tsv")
sysSD1 <- import("results/res_fastlm.sysSD_multiBMI_model.tsv")
diaSD1<- import("results/res_fastlm.diaSD_multiBMI_model.tsv")
pc=import("raw/SCAPIS-DATA-PETITION-507-20230214_principal_components.csv")

pc=pc[,c("Subject","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
names(pc)[1]="subject_id"


data=merge(data,pc,by="subject_id",all.x=T)

colnames(info)[1:3] <- c("id", "name", "level")

sbp1=sbp1[which(sbp1$p.value.infl<0.05 & sign(sbp1$estimate)==sign(sbp1$estimate.infl)),]
dbp1=dbp1[which(dbp1$p.value.infl<0.05 & sign(dbp1$estimate)==sign(dbp1$estimate.infl)),]
sysSD1=sysSD1[which(sysSD1$p.value.infl<0.05 & sign(sysSD1$estimate)==sign(sysSD1$estimate.infl)),]
diaSD1=diaSD1[which(diaSD1$p.value.infl<0.05 & sign(diaSD1$estimate)==sign(diaSD1$estimate.infl)),]

sbp1=sbp1$var.x[which(sbp1$q.value<0.05)]
dbp1 =dbp1$var.x[which(dbp1$q.value<0.05)]
sysSD1=sysSD1$var.x[which(sysSD1$q.value<0.05)]
diaSD1=diaSD1$var.x[which(diaSD1$q.value<0.05)]

# linear regression function

lm.fun <- function(y, x, z,x.var,y.var) {

  tryCatch({

    y <- y
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    prev= (table(data[,"x"]>min(data[,"x"]))["TRUE"]*100)/(length(data[,"x"]))
    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data), message = NA)

  }, warning = function(w) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))

  }, error = function(e) {

    data.frame(var.x=x.var,var.y=y.var,prev=prev,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))

  })

}


#model

mgs <- data[, grep("hMGS.",names(data),value=T)]
cov <- data[, c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI")]
# antibiotics
yi="sbp_m_all"

antib1 <- bplapply(sbp1, function(x) lm.fun(data[which(data$antibio == 0),yi], mgs[which(data$antibio == 0), x], cov[which(data$antibio == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
antib1 <- do.call(rbind, antib1)
names(antib1)[1:2]=c("var.x","var.y")

yi="dbp_m_all"

antib2 <- bplapply(dbp1, function(x) lm.fun(data[which(data$antibio == 0),yi], mgs[which(data$antibio == 0), x], cov[which(data$antibio == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
antib2 <- do.call(rbind, antib2)
names(antib2)[1:2]=c("var.x","var.y")
yi="syssd_all"

antib3 <- bplapply(sysSD1, function(x) lm.fun(data[which(data$antibio == 0),yi], mgs[which(data$antibio == 0), x], cov[which(data$antibio == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
antib3 <- do.call(rbind, antib3)
names(antib3)[1:2]=c("var.x","var.y")
yi="diasd_all"

#antib4 <- bplapply(diaSD1, function(x) lm.fun(data[which(data$antibio == 0),yi], mgs[which(data$antibio == 0), x], cov[which(data$antibio == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
#antib4 <- do.call(rbind, antib4)
#names(antib4)[1:2]=c("var.x","var.y")
antib=rbind(antib1,antib2,antib3)#,antib4)

antib$q.value <- p.adjust(antib$p.value, method = "BH", n = sum(!is.na(antib$p.value)))
colnames(antib)[3:ncol(antib)] <- paste0("antib_", colnames(antib)[3:ncol(antib)])

# ppi

yi="sbp_m_all"

ppi1<- bplapply(sbp1, function(x) lm.fun(data[which(data$ppi == 0),yi], mgs[which(data$ppi == 0), x], cov[which(data$ppi == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
ppi1 <- do.call(rbind, ppi1)
names(ppi1)[1:2]=c("var.x","var.y")
yi="dbp_m_all"

ppi2<- bplapply(dbp1, function(x) lm.fun(data[which(data$ppi == 0),yi], mgs[which(data$ppi == 0), x], cov[which(data$ppi == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
ppi2 <-do.call(rbind, ppi2)
names(ppi2)[1:2]=c("var.x","var.y")
yi="syssd_all"

ppi3<- bplapply(sysSD1, function(x) lm.fun(data[which(data$ppi == 0),yi], mgs[which(data$ppi == 0), x], cov[which(data$ppi == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
ppi3 <- do.call(rbind, ppi3)
names(ppi3)[1:2]=c("var.x","var.y")
yi="diasd_all"

#ppi4<- bplapply(diaSD1, function(x) lm.fun(data[which(data$ppi == 0),yi], mgs[which(data$ppi == 0), x], cov[which(data$ppi == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
#ppi4 <- do.call(rbind, ppi4)
#names(ppi4)[1:2]=c("var.x","var.y")
ppi=rbind(ppi1,ppi2,ppi3)#,ppi4)

ppi$q.value <- p.adjust(ppi$p.value, method = "BH", n = sum(!is.na(ppi$p.value)))
colnames(ppi)[3:ncol(ppi)] <- paste0("ppi_", colnames(ppi)[3:ncol(ppi)])

# crohn
yi="sbp_m_all"

crohn1 <-bplapply(sbp1, function(x) lm.fun(data[which(data$UC == 0),yi], mgs[which(data$UC == 0), x], cov[which(data$UC == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
crohn1 <- do.call(rbind, crohn1)
names(crohn1)[1:2]=c("var.x","var.y")

yi="dbp_m_all"
crohn2 <-bplapply(dbp1, function(x) lm.fun(data[which(data$UC == 0),yi], mgs[which(data$UC == 0), x], cov[which(data$UC == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
crohn2 <- do.call(rbind, crohn2)
names(crohn2)[1:2]=c("var.x","var.y")

yi="syssd_all"
crohn3 <-bplapply(sysSD1, function(x) lm.fun(data[which(data$UC == 0),yi], mgs[which(data$UC == 0), x], cov[which(data$UC == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
crohn3 <- do.call(rbind, crohn3)
names(crohn3)[1:2]=c("var.x","var.y")

#yi="diasd_all"
#crohn4 <-bplapply(diaSD1, function(x) lm.fun(data[which(data$UC == 0),yi], mgs[which(data$UC == 0), x], cov[which(data$UC == 0), ],x,yi), BPPARAM = MulticoreParam(cores))
#crohn4 <- do.call(rbind, crohn4)
#names(crohn4)[1:2]=c("var.x","var.y")

crohn=rbind(crohn1,crohn2,crohn3)#,crohn4)

crohn$q.value <- p.adjust(crohn$p.value, method = "BH", n = sum(!is.na(crohn$p.value)))
colnames(crohn)[3:ncol(crohn)] <- paste0("crohn_", colnames(crohn)[3:ncol(crohn)])

# PCs

cov <- data[, c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

yi="sbp_m_all"

pcs1<-bplapply(sbp1, function(x) lm.fun(data[,yi], mgs[, x], cov,x,yi), BPPARAM = MulticoreParam(cores))
pcs1 <- do.call(rbind, pcs1)
names(pcs1)[1:2]=c("var.x","var.y")

yi="dbp_m_all"
pcs2<-bplapply(dbp1, function(x) lm.fun(data[,yi], mgs[, x], cov,x,yi), BPPARAM = MulticoreParam(cores))
pcs2 <- do.call(rbind, pcs2)
names(pcs2)[1:2]=c("var.x","var.y")

yi="syssd_all"
pcs3<-bplapply(sysSD1, function(x) lm.fun(data[,yi], mgs[, x], cov,x,yi), BPPARAM = MulticoreParam(cores))
pcs3 <- do.call(rbind, pcs3)
names(pcs3)[1:2]=c("var.x","var.y")


pcs=rbind(pcs1,pcs2,pcs3)

pcs$q.value <- p.adjust(pcs$p.value, method = "BH", n = sum(!is.na(pcs$p.value)))
colnames(pcs)[3:ncol(pcs)] <- paste0("pcs_", colnames(pcs)[3:ncol(pcs)])

# make and export table
res=merge(antib,crohn,by=c("var.x","var.y"))
res=merge(res,ppi,by=c("var.x","var.y"))
res=merge(res,pcs,by=c("var.x","var.y"))

res$id=res$var.x


res=merge(info[,c("name","id","level")],res,by="id")
export(res, "results/sensitivity.tsv")

sessionInfo()


