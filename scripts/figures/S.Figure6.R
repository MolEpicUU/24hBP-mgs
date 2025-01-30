rm(list=ls())

library(rio)
library(ggplot2)
library(reshape2)
library(ggforce)
library(ggtext)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
sens <- import("results/sensitivity.tsv")

sbp=import("results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp=import("results/res_fastlm.24DBP_multiBMI_model.tsv")
sbp_sd=import("results/res_fastlm.sysSD_multiBMI_model.tsv")
dbp_sd=import("results/res_fastlm.diaSD_multiBMI_model.tsv")

main <- data.frame(rbind(sbp,dbp,sbp_sd,dbp_sd),stringsAsFactors=F)
main$id=main$var.x

main=main[which(main$p.value.infl<0.05 & sign(main$estimate)==sign(main$estimate.infl)),]
main=main[which(main$q.value<0.05),]


data <- lapply(unique(sens$id), function(id) {
  sens <- sens[which(sens$id == id), ]
  main <- main[which(main$id == id), ]
  estimates <- sens[, which(grepl("estimate", colnames(sens)))]
  #type <- gsub("_estimate", "", colnames(sens[, which(grepl("estimate", colnames(sens)))]))
  n <- sens[, which(grepl("_n", colnames(sens)))]
  estimate <- main[, "estimate"]
  p <- sens[, which(grepl("p.value", colnames(sens)))]
  trait<-sens[,"var.y"]
  diff <- (estimates - estimate) / estimate
  mgs <- paste0(sens$name, " ", sens$id, "\n", "Prevalence = ", round(main$prev , 1), "%")
  prev <- round(main$prev , 1)
  data.frame(mgs, trait,prev, diff, p, n)
})

data <- do.call(rbind, data)

x=melt(data,id.vars=c("mgs","trait","prev"))

x1=cbind(
x[which(x$variable%in%"antib_n"==T),c("mgs","trait","prev","value")],
x[which(x$variable%in%"antib_estimate"==T),"value"],
x[which(x$variable%in%"antib_p.value"==T),"value"])
x1$type="Antib"

x2=cbind(
x[which(x$variable%in%"crohn_n"==T),c("mgs","trait","prev","value")],
x[which(x$variable%in%"crohn_estimate"==T),"value"],
x[which(x$variable%in%"crohn_p.value"==T),"value"])
x2$type="IBD"


x3=cbind(
x[which(x$variable%in%"pcs_n"==T),c("mgs","trait","prev","value")],
x[which(x$variable%in%"pcs_estimate"==T),"value"],
x[which(x$variable%in%"pcs_p.value"==T),"value"])
x3$type="PCs"

x4=cbind(
x[which(x$variable%in%"ppi_n"==T),c("mgs","trait","prev","value")],
x[which(x$variable%in%"ppi_estimate"==T),"value"],
x[which(x$variable%in%"ppi_p.value"==T),"value"])
x4$type="PPI"

names(x1)=names(x2)=names(x3)=names(x4)=c("mgs","trait","prev","n","diff","p","type")

data=data.frame(rbind(x1,x2,x3,x4),stringsAsFactors=F)

data$type <- factor(data$type, levels = unique(c(grep("Antib", data$type, value = T), grep("IBD", data$type, value = T), grep("PPI", data$type, value = T),grep("PCs", data$type, value = T))))
data$p <- ifelse(data$p < 0.05, "*P* value < 0.05", "*P* value >= 0.05")

data$trait=ifelse(data$trait=="sbp_m_all","24-hour SBP",
                 ifelse(data$trait=="dbp_m_all","24-hour DBP",
                        ifelse(data$trait=="syssd_all","Variability of 24-hour SBP",
                               ifelse(data$trait=="diasd_all","Variability of 24-hour DBP",
                                      ifelse(data$trait=="bbps","Office SBP","Office DBP")))))



p1=ggplot(data, aes(x = type, y = diff, color = p)) + geom_hline(yintercept = 0, size = 0.2) + geom_point(size = 2) + facet_wrap_paginate(~ trait+mgs, nrow=3,ncol = 3, scales = "free",drop=T) + lims(y = c(-0.6, 0.6)) + theme_classic() + theme(line = element_line(size = 0.2), text = element_text(size = 9), strip.clip = "off",strip.placement="outside",strip.background = element_blank(), panel.spacing.x = unit(2, "lines"), panel.spacing.y = unit(4, "lines"),legend.text = element_markdown()) + scale_color_manual(values = safe_colorblind_palette[c(8, 10)]) + labs(x = NULL, y = "(sensitivity analysis estimate - model 2 estimate) / model 2 estimate", color = NULL)

pages=n_pages(p1)

pdf("results/figures/S.figure6.pdf", width = 8.3, height = 11.7)
for(i in 1:pages){
print(i)
 print(ggplot(data, aes(x = type, y = diff, color = p)) + geom_hline(yintercept = 0, size = 0.2) + geom_point(size = 2) + facet_wrap_paginate(~ trait+mgs, nrow=3,ncol = 3, scales = "free",drop=T,page=i) + lims(y = c(-0.6, 0.6)) + theme_classic() + theme(line = element_line(size = 0.2), text = element_text(size = 9), strip.clip = "off",strip.placement="outside",strip.background = element_blank(), panel.spacing.x = unit(2, "lines"), panel.spacing.y = unit(4, "lines"),legend.text = element_markdown()) + scale_color_manual(values = safe_colorblind_palette[c(8, 10)]) + labs(x = NULL, y = "(sensitivity analysis estimate - model 2 estimate) / model 2 estimate", color = NULL))
}
dev.off()

