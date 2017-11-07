library(survival)
library(risksetROC)
# Read input data
ads1 <- read.csv('input/170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('input/170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')

# Make ADS
max_therapy_code <- tapply(ads2$REAL_THERAPY_CD, ads2$ID, max)   # extract max therapy code for each patient ID
df1 <- data.frame(ID=row.names(max_therapy_code), max_therapy_code)   # make dataframe only with max therapy code
ads3 <- merge(ads2, df1, by = "ID")
ads4 <- ads3[ads3$REAL_THERAPY_CD == ads3$max_therapy_code, ]   # extract rows only with max therapy code
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
                     ifelse((ads6$CYTO_7 == 1 | ads6$CYTO_5Q == 1 | ads6$CYTO_T1621 == 1 | ads6$CYTO_PH1 == 1 |
                             ads6$FLT3_ITD1 == 2), 3, 2))
ads7 <- ads6[ads6$disease_risk == 1 & !is.na(ads6$ARDI2),]   # extract rows only LR


# ROC analysis, to see how well the marker predicts 3-year survival
auc <- NULL
main_title <- NULL
survival_time <- ads6$EFS_DAY
survival_status <- ads6$EFS_FLG
survival_time1 <- ifelse((ads6$EFS_DAY == 0), 1, ads6$EFS_DAY)
survival_time2 <- ads6$OS_DAY
survival_status2 <- ads6$DETH_FLG

main_title[1] <- "ROC1: Endpoint=EfS, use ARDI, AGE2c, disease_risk"
fit1 <- coxph(Surv(survival_time1, survival_status) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit1)
eta1 <- fit1$linear.predictor
roc1 <- risksetROC(Stime=survival_time1, status=survival_status, marker=eta1, predict.time=1095,
                       method="Cox", main=main_title[1], lty=2, col="red")
out <- CoxWeights(marker=eta1, Stime=survival_time1, status=survival_status, predict.time=1095)
auc[1] <- out$AUC

tmax=1095
AUC.CC1=risksetAUC(Stime=survival_time1, status=survival_status, marker=eta1, tmax=tmax,
                  method="Cox",  lty=2, col="red")

main_title[2] <- "ROC2: Endpoint=OS, use ARDI, AGE2c, disease_risk"
fit2 <- coxph(Surv(survival_time2, survival_status2) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit2)
eta2 <- fit2$linear.predictor
roc2 <- risksetROC(Stime=survival_time2, status=survival_status2, marker=eta2, predict.time=1095,
                       method="Cox", main=main_title[2], lty=2, col="red")
out <- CoxWeights(marker=eta2, Stime=survival_time2, status=survival_status2, predict.time=1095)
auc[2] <- out$AUC

tmax=1095
AUC.CC2=risksetAUC(Stime=survival_time2, status=survival_status2, marker=eta2, tmax=tmax,
                   method="Cox",  lty=2, col="green")


main_title[3] <- "ROC3: Endpoint=OS, use RDI_ANT, AGE2c, disease_risk"
fit3 <- coxph(Surv(survival_time2, survival_status2) ~ RDI_ANT+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit3)
eta3 <- fit3$linear.predictor
roc3 <- risksetROC(Stime=survival_time2, status=survival_status2, marker=eta3, predict.time=1095,
                       method="Cox", main=main_title[3], lty=2, col="red")
out <- CoxWeights(marker=eta3, Stime=survival_time2, status=survival_status2, predict.time=1095)
auc[3] <- out$AUC

tmax=1095
AUC.CC3=risksetAUC(Stime=survival_time2, status=survival_status2, marker=eta3, tmax=tmax,
                   method="Cox",  lty=2, col="blue")

main_title[4] <- "ROC4: Endpoint=OS, use RDI_ARAC, AGE2c, disease_risk"
fit4 <- coxph(Surv(survival_time2, survival_status2) ~ RDI_ARAC+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit4)
eta4 <- fit4$linear.predictor
roc4 <- risksetROC(Stime=survival_time2, status=survival_status2, marker=eta4, predict.time=1095,
                       method="Cox", main=main_title[4], lty=2, col="red")
out <- CoxWeights(marker=eta4, Stime=survival_time2, status=survival_status2, predict.time=1095)
auc[4] <- out$AUC

tmax=1095
AUC.CC4=risksetAUC(Stime=survival_time2, status=survival_status2, marker=eta4, tmax=tmax,
                      method="Cox",  lty=2, col="pink")


main_title[5] <- "ROC5: Endpoint=OS, use RDI_VP16, AGE2c, disease_risk"
fit5 <- coxph(Surv(survival_time2, survival_status2) ~ RDI_VP16+AGE2C+disease_risk, data=ads6, na.action=na.omit)
summary(fit5)
eta5 <- fit5$linear.predictor
roc5 <- risksetROC(Stime=survival_time2, status=survival_status2, marker=eta5, predict.time=1095,
                       method="Cox", main=main_title[5], lty=2, col="red")
out <- CoxWeights(marker=eta5, Stime=survival_time2, status=survival_status2, predict.time=1095)
auc[5] <- out$AUC

AUC.CC5=risksetAUC(Stime=survival_time2, status=survival_status2, marker=eta5, tmax=tmax,
                   method="Cox",  lty=2, col="purple")

main_title[6] <- "ROC6: Endpoint=EfS, use ARDI2 and ads7"
fit6 <- coxph(Surv(survival_time1, survival_status) ~ ARDI, data=ads7, na.action=na.omit)
summary(fit6)
eta6 <- fit6$linear.predictor
roc6 <- risksetROC(Stime=survival_time1, status=survival_status, marker=eta6, predict.time=300,
                   method="Cox", main=main_title[6], lty=2, col="red")
out <- CoxWeights(marker=eta6, Stime=survival_time1, status=survival_status, predict.time=300)
auc[6] <- out$AUC

tmax=300
AUC.CC6=risksetAUC(Stime=survival_time1, status=survival_status, marker=eta6, tmax=tmax,
                   method="Cox",  lty=2, col="red")


