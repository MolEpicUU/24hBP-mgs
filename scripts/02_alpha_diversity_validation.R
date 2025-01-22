rm(list=ls())

print(Sys.time())

#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)

##Folder to save the outputs
output.folder="results/"

# load data
dades=import("0_data/dades.pheno_malmo_2575_uppsala_195_validation.csv")

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
yi="bbps"
#covariates
covari2=c("agev1","gender","q005a","siteid_plate","cur_smoke", "Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa")
covari=c("agev1","gender","q005a","siteid_plate","cur_smoke", "Fibrer","Energi_kcal","diab_treat","HC_treat","sodium_kawa","BMI")
x="shannon"

#run the models

res.lm1 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm1$model="BMI"
res.lm2 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm2$model="multiple"

yi="bbpd"

#run the models

res.lm3 <-fastlm.fun(dades[,x],dades[,yi],dades[,covari],x,yi)
res.lm3$model="BMI"
res.lm4 <- fastlm.fun(dades[,x],dades[,yi],dades[,covari2],x,yi)
res.lm4$model="multiple"

alph.lm=rbind(res.lm1,res.lm2,res.lm3,res.lm4)

export(alph.lm, file = paste(output.folder,"alpha.lm_validation.csv",sep=""), format = "csv")
print(Sys.time())

