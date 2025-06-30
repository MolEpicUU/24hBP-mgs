print(Sys.time())

rm(list=ls())

#load libraries
library(rio)
library(ggplot2)
library(reshape2)
library(compareGroups)
library(ggpubr)
library(Hmisc)
library(vegan)

output.folder="./0_data/"

#load data information
#tax.info=import("raw/HG3.A.7_tax.xlsx")
dades<- import("raw/pheno_MGS_shannon_bray_curtis_MGP_4839_upp_4980_malmo.tsv")
dades.mgs <- import("raw/champ_clr_transformed_combined_data.csv")

dades <- dades[,!names(dades)%in%grep("____",names(dades),value=T)]

dades.bp=import("raw/phenotypes_24h_BP.csv")
dades.con_id=import("raw/id_conversion.txt")
ds=import("raw/scapis_metagenomics_mgs_ds_relative_abundances_v1.0.tsv")
mgs=import("raw/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv")
htn=import("raw/Hypertension_med_12m.tsv")

dades.meta=import("raw/scapis_metabolon_measurements_batchnorm_v2.1.tsv")

drugs <- import("raw/SCAPIS-REGISTERSU-1715-LMED-20220511-T_R_LMED_27947_2021.txt")
sel_participants<-import("raw/selected_participants_id_sample_id.tsv")
dades.diet=import("./raw/AdSugMacro_ScapisMöUppJanuary2021_UE220907_tSergiOct2022_.sav")

sodium.upp=import("./raw/sodium_malmo_uppsala.csv")
sodium.mal=import("./raw/urine_sodium_scapis_lund.csv")

gmm_up=import("raw/upugut03.GMMComp.percent.xlsx")
gmm_m=import("raw/lungut03.GMMComp.percent.xlsx")

names(gmm_up)=gsub("b","",names(gmm_up))
names(gmm_m)=gsub("_b","",names(gmm_m))
names(gmm_m)=gsub("b","",names(gmm_m))
names(gmm_m)=gsub("c","",names(gmm_m))

gmm=merge(gmm_up[which(names(gmm_up)%in%"annotations"==F)],gmm_m[which(names(gmm_m)%in%"annotations"==F)],by="Module")
rownames(gmm)=gmm$Module
gmm=gmm[,grep("gut_",names(gmm),value=T)]
gmm=data.frame(t(gmm),stringsAsFactors=F)
gmm$sample.id=rownames(gmm)
dades=merge(dades,gmm,by="sample.id",all.x=T)

#load a few phenotypes
dades$siteid_plate=paste(dades$siteid,"_",dades$plate,sep="")
dades=dades[,c( "id","day1visitd","siteid","plate","siteid_plate","gender","smokestatus","q005a","bbps","bbpd","agev1","height","weight","bmi","diabd","q030jx", "q030j", "q030kx" , "q030k",  "q030wx", "q030lx", "q030l", "q134")]#,grep("HG3A.",names(dades),value=T),grep("MF",names(dades),value=T))]

#load BP (N=4480)
dades.bp$ABPM="YES"
dades.bp=dades.bp[,c("patno","subject_id", "site","ABPM" , "sbp_m_all", "dbp_m_all", "syssd_all", "diasd_all", "SBP_Mean", "DBP_Mean")]#, "HTN","Height","Weight", "age" ,"Gender" ,"BMI",  "Waist", "Diabetes", "egfr" ,"cur_smoke","CreatinineFormattedResult" ,"area_ses")]

# Load convert ID (N=4800)
dades.bp1 <- merge(dades.bp, dades.con_id, by = "subject_id")

# Merge BP to phenotype (N=4480)
dades.bp1$scapis_rid <- ifelse(dades.bp1$site=="Lund", dades.bp1$subject_id, dades.bp1$export_id)
# merge with microbiome data
dades2=merge(dades.bp1,dades,by.x="scapis_rid",by.y="id",all.y=T)
dades2[which(dades2$ABPM%in%"YES"==F),"ABPM"]="NO"

dades.con_id$subject=ifelse(dades.con_id$subject_id%in%grep("5-",dades.con_id$subject_id,value=T),dades.con_id$export_id,dades.con_id$subject_id)

dades2=merge(dades2[,which(names(dades2)%in%c("subject_id","export_id")==F)],dades.con_id,by.x="scapis_rid",by.y="subject")

#load Metabolomics data(for PPI drug)
dades.meta=dades.meta[,c("scapis_id","met_100002725","met_100002808")]


dades.meta$ppi=ifelse(dades.meta$met_100002725>min(dades.meta$met_100002725,na.rm=T) |
                        dades.meta$met_100002808>min(dades.meta$met_100002808,na.rm=T),1,0)

dades2=merge(dades2,dades.meta,by.x="subject_id",by.y="scapis_id", all.x = TRUE)

# load antibiotics data(for sensitivity test)

narrow <- drugs[which(drugs$ATC %in% c("J01CE02","J01CF05","J01EA01","J01FA01","J01FA06","J01FA09","J01FA10","J01FF01","J01XC01","J01XE01")), ]
names(narrow)=tolower(names(narrow))
narrow$visit <- as.Date(dades2$day1visitd, "%d-%B-%Y")[match(narrow$id, dades2$subject_id)]
narrow$edatum <- as.Date(narrow$edatum)
narrow <- narrow[which(narrow$visit - narrow$edatum < 182 & narrow$visit - narrow$edatum >= 0), ]
dades2$narrow <- ifelse(dades2$subject_id %in% narrow$id, "Yes", "No")
broad <- drugs[which(drugs$ATC %in% c("J01AA02","J01AA04","J01AA06","J01AA07","J01CA04","J01CA08","J01CR02","J01DB05","J01DD14","J01EE01","J01MA02","J01MA06","J01MA12","J01MA14","J01XX05","J01XX08")), ]
names(broad)=tolower(names(broad))
broad$visit <- as.Date(dades2$day1visitd, "%d-%B-%Y")[match(broad$id, dades2$subject_id)]
broad$edatum <- as.Date(broad$edatum)
broad <- broad[which(broad$visit - broad$edatum < 182 & broad$visit - broad$edatum >= 0), ]
dades2$broad <- ifelse(dades2$subject_id %in% broad$id, "Yes", "No")
dades2$antibio=ifelse(dades2$narrow=="Yes" | dades2$broad=="Yes",1,0)
#table(dades2$antibio)

#load diet data(fiber, energy, alcohol)
dades.diet=dades.diet[,c("Subject","Enerkcal","Fibrer")]#,"Alkohol")]

dades2=merge(dades2,dades.diet,by.x="subject_id",by.y="Subject", all.x = TRUE)
dades2=merge(dades2,htn,by.x="subject_id",by.y="ID")

#shannon
ds$shannon <- diversity(ds[,2:ncol(ds)],index="shannon")



#select the variables important for the analysis
var=c("scapis_rid","subject_id", "ABPM","sbp_m_all", "dbp_m_all", "syssd_all", "diasd_all", "bbps", "bbpd", "SBP_Mean", "DBP_Mean","agev1","height","weight","bmi","diabd","smokestatus",# "HTN","Weight","Height","BMI", "Waist", "Diabetes", "egfr" ,"area_ses",
      #"cur_smoke","CreatinineFormattedResult", "age","Gender","HTN",
      "gender","site", "q005a",  "ppi","siteid_plate", "hypermed_register","q030jx", "q030j", "q030kx" , "q030k",  "q030wx", "q030lx", "q030l", "q134","height", "weight", "narrow","broad","Fibrer","Enerkcal","antibio" )#"shannon",grep("___",names(dades2),value=T),grep("MF",names(dades2),value=T))

dades.pheno=dades2[,var]
dades.pheno$age=dades.pheno$agev1
dades.pheno$BMI=dades.pheno$bmi
dades.pheno$Gender=dades.pheno$gender
dades.pheno$Energi_kcal=dades.pheno$Enerkcal

#transform the variable in a factor
dades.pheno$smokestatus=as.factor(dades.pheno$smokestatus)
dades.pheno$cur_smoke=ifelse(dades.pheno$smokestatus==3,1,
                                 ifelse(dades.pheno$smokestatus==2,0,
                                        ifelse(dades.pheno$smokestatus==1,0,
                                               ifelse(dades.pheno$smokestatus==9,NA,1))))
dades.pheno$cur_smoke=as.factor(dades.pheno$cur_smoke)


#Define HBP treatment, high cholesterol, high cholesterol treatment, diabetes and diabetes treatment variables
dades.pheno$HBP_treat=ifelse(dades.pheno$hypermed_register=="no",0,1)

dades.pheno$HC=ifelse(dades.pheno$q030kx==2,NA,
                      ifelse(dades.pheno$q030kx==0,0,1))

dades.pheno$HC_treat=ifelse(dades.pheno$q030kx==0,0,
                            ifelse(dades.pheno$q030k==0,0,
                                   ifelse(dades.pheno$q030k==2,NA,1)))

dades.pheno$UC=ifelse(dades.pheno$q030wx==2,NA,
                      ifelse(dades.pheno$q030wx==0,0,1))


dades.pheno$diab=ifelse(dades.pheno$q030lx==2,NA,
                        ifelse(dades.pheno$q030lx==0,0,1))
dades.pheno$diab_treat=ifelse(dades.pheno$q030lx==0,0,
                              ifelse(dades.pheno$q030l==0,0,
                                     ifelse(dades.pheno$q030l==2,NA,1)))


dades.pheno$PE=ifelse(is.na(dades.pheno$q134),0,dades.pheno$q134)

## Checking outcome variable across site
dades.pheno$siteid=as.factor(dades.pheno$siteid)
xdata= dades.pheno[,c("siteid","dbp_m_all","sbp_m_all")]
xdata=melt(xdata)
xdata$fill=paste(xdata$siteid,xdata$variable,sep="_")

# Handle the missing values in SES
# dades.pheno$ses= ifelse(dades.pheno$area_ses=="Highest", "Highest",
#                         ifelse(dades.pheno$area_ses=="Middle", "Middle",
#                                ifelse(dades.pheno$area_ses=="Lowest", "Lowest", NA)))
#

dades.pheno$siteid=substring(as.character(dades.pheno$siteid),1,1)
dades.pheno$siteid=ifelse(dades.pheno$siteid==2,"Malmö","Uppsala")
dades.pheno$gender=ifelse(dades.pheno$gender==2,"FEMALE","MALE")

#descriptive table
#transform into factor the categorical variables
dades.pheno$siteid_plate=as.factor(dades.pheno$siteid_plate)
dades.pheno$HC=as.factor(dades.pheno$HC)
dades.pheno$diab_treat=as.factor(dades.pheno$diab_treat)
dades.pheno$HC_treat=as.factor(dades.pheno$HC_treat)
dades.pheno$Gender=as.factor(dades.pheno$Gender)
# dades.pheno$HTN=as.factor(dades.pheno$HTN)
dades.pheno$UC=as.factor(dades.pheno$UC)
dades.pheno$ppi=as.factor(dades.pheno$ppi)
dades.pheno$PE=as.factor(dades.pheno$PE)
# dades.pheno$ses=as.factor(dades.pheno$ses)
dades.pheno$antibio=as.factor(dades.pheno$antibio)
dades.pheno$diab=as.factor(dades.pheno$diab)

#sodium
names(sodium.upp)=c("id","sodium","crea","subject_id")
names(sodium.mal)=c("id","sodium","crea")
sodium=data.frame(rbind(sodium.upp[,c("id","sodium","crea")],sodium.mal),stringsAsFactors=F)
dades.pheno<- merge(dades.pheno,sodium, by.x="scapis_rid",by.y="id",all.x=T)
dades.pheno$sodium=ifelse(dades.pheno$sodium=="f?rolyckat"|dades.pheno$sodium<0,NA,ifelse(dades.pheno$sodium=="<3","1.5",ifelse(dades.pheno$sodium=="<4","2",ifelse(dades.pheno$sodium=="<5","2.5",dades.pheno$sodium))))
dades.pheno$crea=ifelse(dades.pheno$crea=="f?rolyckat" | dades.pheno$crea=="f\x9arolyckat",NA,dades.pheno$crea)
dades.pheno$sodium=as.numeric(as.character(dades.pheno$sodium))
dades.pheno$crea=as.numeric(as.character(dades.pheno$crea))
dades.pheno$crea=ifelse(dades.pheno$crea==0,NA,dades.pheno$crea)
dades.pheno$crea=dades.pheno$crea*11.312

#Kawasaki formula
dades.pheno$pre_crea <- ifelse(dades.pheno$gender=="FEMALE",(-4.72*dades.pheno$age)+(8.58*dades.pheno$weight)+(5.09*dades.pheno$height)-74.5, (-12.63*dades.pheno$age)+(15.12*dades.pheno$weight)+(7.39*dades.pheno$height)-79.9 )
dades.pheno$sodium_kawa_original = dades.pheno$sodium_kawa = 23*(16.3*sqrt((dades.pheno$sodium/(dades.pheno$crea*10))*dades.pheno$pre_crea))
dades.pheno$sodium_kawa=log(dades.pheno$sodium_kawa)

dades.pheno=dades.pheno[which(dades.pheno$HBP_treat%in%0==T),]
dades.pheno=dades.pheno[which(is.na(dades.pheno$q005a)==F),]

dades.abpm=dades.pheno[which(dades.pheno$ABPM%in%"YES"==T),]
dades.non=dades.pheno[which(dades.pheno$ABPM%in%"YES"==F),]
dades.non=dades.non[!(is.na(dades.non$bbps) & is.na(dades.non$bbpd)),]

covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")

dades.abpm=dades.abpm[complete.cases(dades.abpm[,covari]),]
dades.non=dades.non[complete.cases(dades.non[,covari]),]


##################################################################
#MGS filter <1% non-zero values (microbiome)
#select only mgs columns
#dades.mgs=dades.abpm[,grep("___",names(dades.abpm),value=T)]
#apply the filter
perc=nrow(dades.abpm)*30/100

mgs.abpm=mgs[which(mgs$scapis_id%in%dades.abpm$subject_id),]
mgs.filter.abpm=apply(mgs.abpm[,2:ncol(mgs.abpm)], 2, function(x) table(x>0.01)["TRUE"]<=perc)
#show how mgs will be removed (will be removed the true mgs)
#table(mgs.filter)
#select the mgs that will be selected
mgs.filter.abpm=names(mgs.filter.abpm[which(mgs.filter.abpm%in%"FALSE")])
#QC checks
#select the selected mgs from dades.mgs
dades.mgs=dades.mgs[,names(dades.mgs)%in%c("subject_id",mgs.filter.abpm)]
#check that any mgs will sum 0
#table(colSums(dades.mgs)==0)


# merge with mgs data

dades.abpm=merge(dades.abpm,dades.mgs,by="subject_id")
dades.non=merge(dades.non,dades.mgs,by="subject_id")

names(ds)[1]="subject_id"
dades.abpm=merge(dades.abpm,ds[,c("subject_id","shannon")],by="subject_id")
dades.non=merge(dades.non,ds[,c("subject_id","shannon")],by="subject_id")

#covari=c("age","gender","q005a","siteid_plate","cur_smoke","Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")

#dades.abpm=dades.abpm[complete.cases(dades.abpm[,covari]),]
#dades.non=dades.non[complete.cases(dades.non[,covari]),]

#check that remaining mgs will have >30 non zero values
#apply the filter
#dades.mgs=dades.abpm[,grep("^hMGS.",names(dades.abpm),value=T)]
#mgs.filter2=apply(dades.mgs[,2:ncol(dades.mgs)], 2, function(x) table(x>0.01)["TRUE"]<perc)
#mgs.filter2=names(mgs.filter2[which(mgs.filter2==TRUE)])


#remove the remianing rare mgs
#try(dades.abpm=dades.abpm[,!names(dades.abpm)%in%mgs.filter2])
#try(dades.non=dades.non[,!names(dades.non)%in%mgs.filter2])
##################################################################



#saving the final phenotype data.frame
write.table(dades.abpm, file=paste(output.folder,'/dades.pheno_malmo_758_uppsala_2937_new.csv',sep=""),
            col.names=T,row.names=F,sep=",")

write.table(dades.non, file=paste(output.folder,'/dades.pheno_malmo_2575_uppsala_195_validation.csv',sep=""),
            col.names=T,row.names=F,sep=",")

print(Sys.time())


