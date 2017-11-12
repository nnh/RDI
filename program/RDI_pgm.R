library(survival)
library(risksetROC)
# Read input data
ads1 <- read.csv('input/170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('input/170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')

# Make ADS
max_therapy_code <- tapply(ads2$REAL_THERAPY_CD, ads2$ID, max)   # extract max therapy code for each patient ID
df1 <- data.frame(ID=row.names(max_therapy_code), max_therapy_code)   # make dataframe only with max therapy code
ads3 <- merge(ads2, df1, by = "ID")
ads6 <- ads3[ads3$REAL_THERAPY_CD == ads3$max_therapy_code & !is.na(ads3$EFF), ]   # only with max therapy code
# ads5 <- ads4[!is.na(ads4$EFS_DAY) & !is.na(ads4$EFS_FLG), ]   # extract rows w/o NA in EFS_DAY and EFS_FLG
# ads6 <- ads4[!is.na(ads4$EFF), ]   # extract rows w/o NA in EFF (efficacy analysis set)

ads6$CYTO_T821[is.na(ads6$CYTO_T821)] <- 0
ads6$CYTO_INV16[is.na(ads6$CYTO_INV16)] <- 0
ads6$FLT3_ITD1[is.na(ads6$FLT3_ITD1)] <- 0
ads6$CYTO_7[is.na(ads6$CYTO_7)] <- 0
ads6$CYTO_5Q[is.na(ads6$CYTO_5Q)] <- 0
ads6$CYTO_T1621[is.na(ads6$CYTO_T1621)] <- 0
ads6$CYTO_PH1[is.na(ads6$CYTO_PH1)] <- 0
ads6$FLT3_ITD1[is.na(ads6$FLT3_ITD1)] <- 0
ads6$EFS_1 <- ifelse((ads6$EFS_DAY == 0), 1, ads6$EFS_DAY)
ads6$disease_risk <- ifelse((ads6$CYTO_T821 == 1 | ads6$CYTO_INV16 == 1) & (ads6$FLT3_ITD1 == 1), 1,
                     ifelse((ads6$CYTO_7 == 1 | ads6$CYTO_5Q == 1 | ads6$CYTO_T1621 == 1 | ads6$CYTO_PH1 == 1 |
                             ads6$FLT3_ITD1 == 2), 3, 2))
ads7 <- ads6[ads6$disease_risk == 1 & !is.na(ads6$ARDI2), ]   # extract rows only LR

# ROC analysis, to see how well the marker predicts 3-year survival
auc <- NULL
main_title <- NULL
predict_time <- NULL
span_ads6_os <- 1.0 * (length(ads6$OS_DAY[ads6$DETH_FLG == 1]) ^ -0.2)
span_ads6_efs <- 1.0 * (length(ads6$EFS_1[ads6$EFS_FLG == 1]) ^ -0.2)
span_ads7_efs <- 1.0 * (length(ads7$EFS_1[ads7$EFS_FLG == 1]) ^ -0.2)

main_title[1] <- "ROC/AUS1: EFS, use ARDI, AGE2c, disease_risk"
predict_time <- 1095
fit1 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta1 <- fit1$linear.predictors
roc1 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta1, predict.time=predict_time,
                   method="LocalCox", span=span_ads6_efs, prop=1.0, main=main_title[1], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc1$AUC,3)), cex=1.5)
auc1 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta1, tmax=predict_time,
                   method="LocalCox", span=span_ads6_efs, main=main_title[1], lty=2, col="red")

main_title[2] <- "ROC/AUS2: OS, use ARDI, AGE2c, disease_risk"
predict_time <- 1095
fit2 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta2 <- fit2$linear.predictors
roc2 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta2, predict.time=predict_time,
                   method="LocalCox", span=span_ads6_os, prop=1.0, main=main_title[2], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc2$AUC,3)), cex=1.5)
auc2 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta2, tmax=predict_time,
                   method="LocalCox", span=span_ads6_os, main=main_title[2], lty=2, col="green")

main_title[3] <- "ROC/AUS3: OS, use RDI_ANT, AGE2c, disease_risk"
predict_time <- 1095
fit3 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_ANT+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta3 <- fit3$linear.predictors
roc3 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta3, predict.time=predict_time,
                   method="LocalCox", span=span_ads6_os, prop=1.0, main=main_title[3], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc3$AUC,3)), cex=1.5)
auc3 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta3, tmax=predict_time,
                   method="LocalCox", span=span_ads6_os, main=main_title[3], lty=2, col="blue")

main_title[4] <- "ROC/AUS4: OS, use RDI_ARAC, AGE2c, disease_risk"
predict_time <- 1095
fit4 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_ARAC+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta4 <- fit4$linear.predictors
roc4 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta4, predict.time=predict_time,
                   method="LocalCox", span=span_ads6_os, prop=1.0, main=main_title[4], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc4$AUC,3)), cex=1.5)
auc4 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta4, tmax=predict_time,
                   method="LocalCox", span=span_ads6_os, main=main_title[4], lty=2, col="pink")

main_title[5] <- "ROC/AUS5: OS, use RDI_VP16, AGE2c, disease_risk"
predict_time <- 1095
fit5 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_VP16+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta5 <- fit5$linear.predictors
roc5 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta5, predict.time=predict_time,
                   method="LocalCox", span=span_ads6_os, prop=1.0, main=main_title[5], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc5$AUC,3)), cex=1.5)
auc5 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta5, tmax=predict_time,
                   method="LocalCox", span=span_ads6_os, main=main_title[5], lty=2, col="purple")

main_title[6] <- "ROC/AUS6: EFS, use ARDI2 and ads7"
predict_time <- 300
fit6 <- coxph(Surv(ads7$EFS_1, ads7$EFS_FLG) ~ ARDI2, data=ads7, na.action=na.omit)
eta6 <- fit6$linear.predictors
roc6 <- risksetROC(Stime=ads7$EFS_1, status=ads7$EFS_FLG, marker=eta6, predict.time=predict_time,
                   method="LocalCox", span=span_ads7_efs, prop=1.0, main=main_title[6], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc6$AUC,3)), cex=1.5)
auc6 <- risksetAUC(Stime=ads7$EFS_1, status=ads7$EFS_FLG, marker=eta6, tmax=predict_time,
                   method="LocalCox", span=span_ads7_efs, main=main_title[6], lty=2, col="red")

main_title[7] <- "ROC/AUC7: EFS, use ARDI"
predict_time <- 1095
fit7 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ ARDI, data=ads6, na.action=na.omit)
eta7 <- fit7$linear.predictors
roc7 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta7, predict.time=predict_time,
                   method="LocalCox", plot=TRUE, span=span_ads6_efs, prop=1.0, main=main_title[7], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc7$AUC,3)), cex=1.5)
auc7 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta7, tmax=predict_time,
                   method="LocalCox", span=span_ads6_efs, main=main_title[7], lty=2, col="red")

main_title[8] <- "ROC/AUC8: OS, use ARDI1"
predict_time <- 100
fit8 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ ARDI1, data=ads6, na.action=na.omit)
eta8 <- fit8$linear.predictors
roc8 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta8, predict.time=predict_time,
                   method="LocalCox", plot=TRUE, span=span_ads6_efs, prop=1.0, main=main_title[8], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc8$AUC,3)), cex=1.5)
auc8 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta8, tmax=predict_time,
                   method="LocalCox", span=span_ads6_os, main=main_title[8], lty=2, col="red")

# cox.zph(fit4)
