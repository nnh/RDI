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

ads6$CYTO_T821[is.na(ads6$CYTO_T821)] <- 0
ads6$CYTO_INV16[is.na(ads6$CYTO_INV16)] <- 0
ads6$FLT3_ITD1[is.na(ads6$FLT3_ITD1)] <- 0
ads6$CYTO_7[is.na(ads6$CYTO_7)] <- 0
ads6$CYTO_5Q[is.na(ads6$CYTO_5Q)] <- 0
ads6$CYTO_T1621[is.na(ads6$CYTO_T1621)] <- 0
ads6$CYTO_PH1[is.na(ads6$CYTO_PH1)] <- 0
ads6$FLT3_ITD1[is.na(ads6$FLT3_ITD1)] <- 0

ads6$disease_risk <- ifelse((ads6$CYTO_T821 == 1 | ads6$CYTO_INV16 == 1) & (ads6$FLT3_ITD1 == 1), 1,
                     ifelse((ads6$CYTO_7 == 1 | ads6$CYTO_5Q == 1 | ads6$CYTO_T1621 == 1 | ads6$CYTO_PH1 == 1 | ads6$FLT3_ITD1 == 2), 3, 2))

# ROC analysis
survival.time <- ads6$EFS_DAY
survival.status <- ads6$EFS_FLG
survival.time1 <- ifelse((ads6$EFS_DAY == 0),1,(ads6$EFS_DAY))
survival.time2 <- ads6$OS_DAY
survival.status2 <- ads6$DETH_FLG


fit0 <- coxph(Surv(survival.time1, survival.status) 
              ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit0)
eta <- fit0$linear.predictor
ROC.CC10 <- risksetROC(Stime=survival.time1, status=survival.status, marker=eta, predict.time=1095, 
                       method="Cox", main="ROC Curve", lty=2, col="red")

AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time1, status=survival.status, predict.time=1095)
                                       # to see how well the marker predicts 3-year survival
AUC <- out$AUC
