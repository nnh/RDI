library(survival)
library(risksetROC)
# Read input data
ads1 <- read.csv('input/170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('input/170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')

# Make ADS
MAX_T <- tapply(ads2$REAL_THERAPY_CD, ads2$ID, max)   # extract max therapy code for each patient ID
max_therapy <- data.frame(MAX_T)   # make dataframe only with max therapy code
max_therapy$ID <- row.names(MAX_T)   # make ID column
ads3 <- merge(ads2, max_therapy, by = "ID")
ads4 <- ads3[ads3$REAL_THERAPY_CD == ads3$MAX_T, ]   # extract rows only with max therapy code
ads5 <- ads4[!is.na(ads4$EFS_DAY) & !is.na(ads4$EFS_FLG), ]   # extract rows w/o NA in EFS_DAY and EFS_FLG
ads6 <- ads4[!is.na(ads4$EFF), ]   # extract rows w/o NA in EFF (efficacy analysis set)
ads7 <- ads6[!is.na(ads6$FLT3_ITD1) & !is.na(ads6$SEX) & !is.na(ads6$AGE10) & !is.na(ads6$CYTO_FU), ]   # extract rows w/o FLT3_ITD1 and SEX and AGE1 and CYTO_FU from EFF pop. and 

# ROC analysis
survival.time <- ads7$EFS_DAY
survival.status <- ads7$EFS_FLG
survival.time2 <- ads7$OS_DAY
survival.status2 <- ads7$DETH_FLG

fit0 <- coxph(Surv(survival.time2, survival.status2) ~ ARDI+FLT3_ITD1+CYTO_FU+RISK, data=ads7, na.action=na.omit)
eta <- fit0$linear.predictors
ROC.CC10 <- risksetROC(Stime=survival.time2, status=survival.status2, marker=eta, predict.time=1095,
                       method="Cox", main="ROC Curve", lty=2, col="red")
AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time2, status=survival.status2, predict.time=1095)
AUC <- out$AUC
