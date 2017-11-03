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

# ROC analysis
survival.time <- ads5$EFS_DAY
survival.status <- ads5$EFS_FLG

fit0 <- coxph(Surv(EFS_DAY, EFS_FLG) ~ ARDI1, data=ads5, na.action=na.omit)
eta <- fit0$linear.predictors
ROC.CC10 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=180,
                       method="Cox", main="ROC Curve", lty=2, col="red")
AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status, predict.time=180)
AUC <- out$AUC
