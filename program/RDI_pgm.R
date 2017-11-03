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
ads5 <- ads4[!is.na(ads4$EFS_DAY) & !is.na(ads4$EFS_FLG), ]   # extract rows without NA in EFS_DAY and EFS_FLG

ads5$CYTO_T821[is.na(ads5$CYTO_T821)] <- 0
ads5$CYTO_INV16[is.na(ads5$CYTO_INV16)] <- 0
ads5$FLT3_ITD1[is.na(ads5$FLT3_ITD1)] <- 0
ads5$CYTO_7[is.na(ads5$CYTO_7)] <- 0
ads5$CYTO_5Q[is.na(ads5$CYTO_5Q)] <- 0
ads5$CYTO_T1621[is.na(ads5$CYTO_T1621)] <- 0
ads5$CYTO_PH1[is.na(ads5$CYTO_PH1)] <- 0
ads5$FLT3_ITD1[is.na(ads5$FLT3_ITD1)] <- 0

ads5$disease_risk <- ifelse((ads5$CYTO_T821 == 1 | ads5$CYTO_INV16 == 1) & (ads5$FLT3_ITD1 == 1), 1,
                     ifelse((ads5$CYTO_7 == 1 | ads5$CYTO_5Q == 1 | ads5$CYTO_T1621 == 1 | ads5$CYTO_PH1 == 1 | ads5$FLT3_ITD1 == 1), 3, 2))

# ROC analysis
survival.time <- ads5$EFS_DAY
survival.status <- ads5$EFS_FLG

fit0 <- coxph(Surv(EFS_DAY, EFS_FLG) ~ ARDI, data=ads5, na.action=na.omit)
eta <- fit0$linear.predictors
ROC.CC10 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=10, method="Cox",
                       main="ROC Curve", lty=2, col="red")
AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status, predict.time=10)
AUC <- out$AUC

# REFERNCE
# library(survival)
# library(risksetROC)
# library(MASS)
# data(VA)
# survival.time <- VA$stime
# survival.status <- VA$status
# score <- VA$Karn
# cell.type <- factor(VA$cell)
# tx <- as.integer(VA$treat == 1)
# age <- VA$age
# survival.status[survival.time > 500 ] <- 0
# survival.time[survival.time > 500 ] <- 500
# fit0 <- coxph(Surv(survival.time, survival.status) ~ score + cell.type + tx + age, na.action=na.omit)
# summary(fit0)
# eta <- fit0$linear.predictor
# ROC.CC30 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=30,
#                        method="Cox", main="ROC Curve", lty=2, col="red")
# AUC <- NULL
# out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status, predict.time=30)
# AUC <- out$AUC   ## to see how well the marker predicts one-month survival

