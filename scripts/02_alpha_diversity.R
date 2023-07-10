rm(list=ls())

print(Sys.time())

#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)

##Folder to save the outputs
output.folder="results/"

# load data
dades=import("0_data/dades.pheno_malmo_802_uppsala_3261_new.csv")

####################################

######################

# regression functions
fastlm.fun <- function(x,y, z,x.var,y.var) {
  
  tryCatch({
    
    y <- y
    x <- x
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(var.x=x.var,var.y=y.var,estimate = coef["x", 1], lower = ci["x", 1], upper = ci["x", 2], se = coef["x", 2], p.value = coef["x", 4], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(var.x=x.var,var.y=y.var,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(var.x=x.var,var.y=y.var,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

#outcome
yi="sbp_m_all"
#covariates
covari=c("age","gender", "q005a","siteid_plate")
covari2=c("age","gender","q005a","siteid_plate","cur_smoke", "Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")
covari3=c("age","gender","q005a","siteid_plate","cur_smoke", "Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI")
x="shannon"

#run the models

res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm1$model="crude"
res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm2$model="multiple"
res.lm3 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm3$model="BMI"

yi="dbp_m_all"

#run the models

res.lm4 <-fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm4$model="crude"
res.lm5 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm5$model="multiple"
res.lm6 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm6$model="BMI"

yi="syssd_all"

#run the models

res.lm7 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm7$model="crude"
res.lm8 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm8$model="multiple"
res.lm9 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm9$model="BMI"

yi="diasd_all"

#run the models

res.lm10 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm10$model="crude"
res.lm11 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm11$model="multiple"
res.lm12 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm12$model="BMI"
############################
##########office BP#########
############################


yi="bbps"
#run the models

res.lm13 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm13$model="crude"
res.lm14 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm14$model="multiple"
res.lm15 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm15$model="BMI"


yi="bbpd"

#run the models

res.lm16 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm16$model="crude"
res.lm17 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm17$model="multiple"
res.lm18 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari3],x,yi)
res.lm18$model="BMI"


alph.lm=rbind(res.lm1,res.lm2,res.lm3,res.lm4,res.lm5,res.lm6,res.lm7,res.lm8,res.lm9,res.lm10,res.lm11,res.lm12,res.lm13,res.lm14,res.lm15,res.lm16,res.lm17,res.lm18)

export(alph.lm, file = paste(output.folder,"alpha.lm.csv",sep=""), format = "csv")
print(Sys.time())

