rm(list=ls())
#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)
library(reshape2)
library(ComplexHeatmap)

#load results
sbp1=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.24SBP_multiBMI_model.tsv")
dbp1=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.24DBP_multiBMI_model.tsv")
sbp_sd1=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.sysSD_multiBMI_model.tsv")
dbp_sd1=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.diaSD_multiBMI_model.tsv")
sbp2=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.24SBP_multi_model.tsv")
dbp2=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.24DBP_multi_model.tsv")
sbp_sd2=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.sysSD_multi_model.tsv")
dbp_sd2=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/res_fastlm.diaSD_multi_model.tsv")

info=import("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/MolEpic/Projects/Ongoing_Projects/Blood_pressure_Microbiota_Emily/revision_1_CHAMP/results/raw/scapis_metagenomics_mgs_annotations_v1.0.tsv")
names(info)[1]="var.x"

sel=unique(c(sbp1[which(sbp1$q.value<0.05),"var.x"],dbp1[which(dbp1$q.value<0.05),"var.x"],sbp_sd1[which(sbp_sd1$q.value<0.05),"var.x"],dbp_sd1[which(dbp_sd1$q.value<0.05),"var.x"],
             sbp2[which(sbp2$q.value<0.05),"var.x"],dbp2[which(dbp2$q.value<0.05),"var.x"],sbp_sd2[which(sbp_sd2$q.value<0.05),"var.x"],dbp_sd2[which(dbp_sd2$q.value<0.05),"var.x"]))

prev=merge(info,sbp2[,c("prev","var.x")],by="var.x")
prev=prev[which(prev$var.x%in%sel),]
prev$prev=round(prev$prev,1)
prev$var.x=paste(prev$maintax," (",prev$var.x,")",sep="")



sbp1=sbp1[which(sbp1$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
sbp1$q.value=ifelse(sbp1$q.value<0.05,"**",ifelse(sbp1$p.value<0.05,"*",""))
dbp1=dbp1[which(dbp1$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
dbp1$q.value=ifelse(dbp1$q.value<0.05,"**",ifelse(dbp1$p.value<0.05,"*",""))
sbp_sd1=sbp_sd1[which(sbp_sd1$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
sbp_sd1$q.value=ifelse(sbp_sd1$q.value<0.05,"**",ifelse(sbp_sd1$p.value<0.05,"*",""))
dbp_sd1=dbp_sd1[which(dbp_sd1$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
dbp_sd1$q.value=ifelse(dbp_sd1$q.value<0.05,"**",ifelse(dbp_sd1$p.value<0.05,"*",""))

sbp2=sbp2[which(sbp2$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
sbp2$q.value=ifelse(sbp2$q.value<0.05,"**",ifelse(sbp2$p.value<0.05,"*",""))
dbp2=dbp2[which(dbp2$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
dbp2$q.value=ifelse(dbp2$q.value<0.05,"**",ifelse(dbp2$p.value<0.05,"*",""))
sbp_sd2=sbp_sd2[which(sbp_sd2$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
sbp_sd2$q.value=ifelse(sbp_sd2$q.value<0.05,"**",ifelse(sbp_sd2$p.value<0.05,"*",""))
dbp_sd2=dbp_sd2[which(dbp_sd2$var.x%in%sel),c("var.x","var.y","prev","estimate","q.value","p.value")]
dbp_sd2$q.value=ifelse(dbp_sd2$q.value<0.05,"**",ifelse(dbp_sd2$p.value<0.05,"*",""))

sbp1=merge(sbp1,info[,c("var.x","maintax")],by="var.x")
dbp1=merge(dbp1,info[,c("var.x","maintax")],by="var.x")
sbp_sd1=merge(sbp_sd1,info[,c("var.x","maintax")],by="var.x")
dbp_sd1=merge(dbp_sd1,info[,c("var.x","maintax")],by="var.x")

sbp2=merge(sbp2,info[,c("var.x","maintax")],by="var.x")
dbp2=merge(dbp2,info[,c("var.x","maintax")],by="var.x")
sbp_sd2=merge(sbp_sd2,info[,c("var.x","maintax")],by="var.x")
dbp_sd2=merge(dbp_sd2,info[,c("var.x","maintax")],by="var.x")

sbp1$var.x=paste(sbp1$maintax," (",sbp1$var.x,")",sep="")
dbp1$var.x=paste(dbp1$maintax," (",dbp1$var.x,")",sep="")
sbp_sd1$var.x=paste(sbp_sd1$maintax," (",sbp_sd1$var.x,")",sep="")
dbp_sd1$var.x=paste(dbp_sd1$maintax," (",dbp_sd1$var.x,")",sep="")

sbp2$var.x=paste(sbp2$maintax," (",sbp2$var.x,")",sep="")
dbp2$var.x=paste(dbp2$maintax," (",dbp2$var.x,")",sep="")
sbp_sd2$var.x=paste(sbp_sd2$maintax," (",sbp_sd2$var.x,")",sep="")
dbp_sd2$var.x=paste(dbp_sd2$maintax," (",dbp_sd2$var.x,")",sep="")


sbp1=na.omit(sbp1)
dbp1=na.omit(dbp1)
sbp_sd1=na.omit(sbp_sd1)
dbp_sd1=na.omit(dbp_sd1)
sbp2=na.omit(sbp2)
dbp2=na.omit(dbp2)
sbp_sd2=na.omit(sbp_sd2)
dbp_sd2=na.omit(dbp_sd2)

res.sbp1 <- dcast(sbp1,var.x~var.y,value.var="estimate")
res.sbp1.q <- dcast(sbp1,var.x~var.y,value.var="q.value")
res.dbp1 <- dcast(dbp1,var.x~var.y,value.var="estimate")
res.dbp1.q <- dcast(dbp1,var.x~var.y,value.var="q.value")
res.sbp_sd1 <- dcast(sbp_sd1,var.x~var.y,value.var="estimate")
res.sbp_sd1.q <- dcast(sbp_sd1,var.x~var.y,value.var="q.value")
res.dbp_sd1 <- dcast(dbp_sd1,var.x~var.y,value.var="estimate")
res.dbp_sd1.q <- dcast(dbp_sd1,var.x~var.y,value.var="q.value")

res.sbp2 <- dcast(sbp2,var.x~var.y,value.var="estimate")
res.sbp2.q <- dcast(sbp2,var.x~var.y,value.var="q.value")
res.dbp2 <- dcast(dbp2,var.x~var.y,value.var="estimate")
res.dbp2.q <- dcast(dbp2,var.x~var.y,value.var="q.value")
res.sbp_sd2 <- dcast(sbp_sd2,var.x~var.y,value.var="estimate")
res.sbp_sd2.q <- dcast(sbp_sd2,var.x~var.y,value.var="q.value")
res.dbp_sd2 <- dcast(dbp_sd2,var.x~var.y,value.var="estimate")
res.dbp_sd2.q <- dcast(dbp_sd2,var.x~var.y,value.var="q.value")

res1=merge(prev[,c("var.x","prev")],res.sbp1,by="var.x",all.y=T)
res1=merge(res1,res.dbp1,by="var.x",all.x=T,all.y=T)
res1=merge(res1,res.sbp_sd1,by="var.x",all.x=T,all.y=T)
res1=merge(res1,res.dbp_sd1,by="var.x",all.x=T,all.y=T)

res1=merge(res1,res.sbp2,by="var.x",all.x=T,all.y=T)
res1=merge(res1,res.dbp2,by="var.x",all.x=T,all.y=T)
res1=merge(res1,res.sbp_sd2,by="var.x",all.x=T,all.y=T)
res1=merge(res1,res.dbp_sd2,by="var.x",all.x=T,all.y=T)

res1.q=merge(res.sbp1.q,res.dbp1.q,by="var.x",all.x=T,all.y=T)
res1.q=merge(res1.q,res.sbp_sd1.q,by="var.x",all.x=T,all.y=T)
res1.q=merge(res1.q,res.dbp_sd1.q,by="var.x",all.x=T,all.y=T)

res1.q=merge(res1.q,res.sbp2.q,by="var.x",all.x=T,all.y=T)
res1.q=merge(res1.q,res.dbp2.q,by="var.x",all.x=T,all.y=T)
res1.q=merge(res1.q,res.sbp_sd2.q,by="var.x",all.x=T,all.y=T)
res1.q=merge(res1.q,res.dbp_sd2.q,by="var.x",all.x=T,all.y=T)

res1[is.na(res1)]=0
res1.q[is.na(res1.q)]=""

# Row annotation 
rowannot=res1$prev
names(rowannot)=res1$var.x

rowannot=rowAnnotation(Prevalence = anno_barplot(rowannot, border = F, axis_param = list(gp=gpar(fontsize=6))), annotation_name_side = "top", 
                       annotation_name_gp = gpar(fontsize=6))

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#6E8B3D", "#6E8B3D","#6E8B3D","#6E8B3D"),
                                        labels = c("SBP", "DBP", "SBP-SD","DBP-SD"),
                                        labels_gp = gpar(col = "white", fontsize = 7)),
                       height = unit(0.3, "cm"))

names(res1)=c("var.x","prev","mod1_sbp","mod1_dbp","mod1_sbp_sd","mod1_dbp_sd","mod2_sbp","mod2_dbp",
              "mod2_sbp_sd","mod2_dbp_sd")

names(res1.q)=c("var.x","mod1_sbp","mod1_dbp","mod1_sbp_sd","mod1_dbp_sd","mod2_sbp","mod2_dbp",
                "mod2_sbp_sd","mod2_dbp_sd")

res1=res1[,c("var.x","prev","mod1_sbp","mod2_sbp","mod1_dbp","mod2_dbp",
             "mod1_sbp_sd","mod2_sbp_sd","mod1_dbp_sd","mod2_dbp_sd")]
res1.q=res1.q[,c("var.x","mod1_sbp","mod2_sbp","mod1_dbp","mod2_dbp",
                 "mod1_sbp_sd","mod2_sbp_sd","mod1_dbp_sd","mod2_dbp_sd")]


matrix.res1=as.matrix(res1[,3:ncol(res1)])
rownames(matrix.res1)=res1$var.x
res1.q=res1.q[,2:ncol(res1.q)]
rownames(res1.q)=res1.q$var.x

split2 = rep(1:4, each = 2)
split2=ifelse(split2==1," ",ifelse(split2==2,"  ",ifelse(split2==3,"   ","    ")))

h1=  Heatmap(matrix.res1,
             column_split = split2,
             top_annotation = ha,
             left_annotation = rowannot,
             column_names_rot =  30,
             column_names_side = "top",
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             column_title = "",
             #name="BPs",
             show_heatmap_legend = FALSE,
             column_labels = c("Model 1", "Model 2", "Model 1", "Model 2", "Model 1", "Model 2", "Model 1", "Model 2"),
             column_names_gp = gpar(fontsize = 8),
             
             #column_title_gp = gpar(fontsize = 7,fontface='bold'),
             #heatmap_legend_param = list(labels_gp = gpar(fontsize=7),
             #                            title_gp = gpar(fontsize=8),
             #                            grid_width = unit(5, "mm")),
             #left_annotation = ha,
             show_row_dend = FALSE,
             row_names_side = "left",
             row_names_gp = gpar(fontsize = 5),
             col = colorRamp2(c(-0.5,0,0.5),c("darkblue","white","darkred")),
             cell_fun = function(j, i, x, y, width, height, fill)
             {
               gb = textGrob("*")
               gb_w = convertWidth(grobWidth(gb), "mm")
               gb_h = convertHeight(grobHeight(gb), "mm")
               grid.text(sprintf("%s", res1.q[i, j]),
                         x, y - gb_h*0.4 + gb_w*0.4, gp = gpar(fontsize = 7))
             },
             row_dend_width = unit(1.3, "cm"),
             clustering_distance_rows = "euclidean",
             #layer_fun = function(j, i, x, y, w, h, fill){
             #   pindex(hm.matrix1,j, i)
             #   grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
             #  }
             width=ncol(matrix.res1)*unit(8,"mm"),
             height=nrow(matrix.res1)*unit(1.5,"mm")
             
)

# pdf("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/Fig1.pdf")
draw(h1)
# dev.off()

