library(rio)
x=import("/proj/nobackup/sens2019512/Projects/24hbp_mgs/0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")
y=import("/proj/nobackup/sens2019512/Projects/24hbp_mgs/0_data/dades.pheno_malmo_2866_uppsala_395_validation.csv")
z=data.frame(rbind(x,y),stringsAsFactors=F)
write.table(z,file="/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1/0_combine_data.tsv",sep="\t",col.names=T,row.names=F)

