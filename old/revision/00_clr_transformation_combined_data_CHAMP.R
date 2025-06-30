# Transformation of the microbiome data for the 24hbp project

library(data.table)

setwd("/proj/nobackup/sens2019512/Projects/24hbp_mgs/revision_1/")

pheno <- fread('0_combine_data.tsv')
champ <- fread('/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv', data.table=F)

  # CLR transformation function:
  clr_fun <- function(spmat){
    
    zeros <- lapply(1:ncol(spmat), function(i) which(spmat[,i]==0))
    
    min_non_zero <- min(spmat[spmat>0])
    
    tempmat <- spmat + min_non_zero

    
    mgsclr <- as.matrix(as.data.frame(compositions::clr(tempmat)))
    
    # Replacing the inicial zeros with the minimal value for every species.
    for(i in 1:ncol(mgsclr)){
      
      min <- min(mgsclr[,i][ -zeros[[i]] ])
      mgsclr[,i][ zeros[[i]]  ] <- min
    }
    data.frame(subject_id = rownames(mgsclr), mgsclr)
    
  }
  
  # Transformation  ---------------------------------------------------------------------

  champ <- champ[champ$scapis_id %in% pheno$subject_id, ]
  
  rownames(champ) <- champ$scapis_id
  
  allzeros    <- names(champ[,-1])[colSums(champ[,-1])==0]
  champ  <- champ[, -which(colnames(champ) %in% allzeros)]
  
  champ  <- as.matrix(champ[,-1])
  
  
  champ_clr <- clr_fun(champ)
  
  message("\nCLR Transformation finished")
  
  
  fwrite(champ_clr, file = 'champ_clr_transformed_combined_data.csv')
  
  
  message("\nEnd")
