library(survival)
library(risksetROC)

OpenOutputdevice <- function(extension, filepath){
  if (extension == "eps") {
    setEPS()
    postscript(filepath)
  } else if (extension == "png") {
    png(filepath)
  } else if (extension == "wmf") {
    win.metafile(filename=filepath)
  }
}

RocAuc2 <- function(extension, predict_time, survival_time, status, event, variable, dataframe){
  fit <- coxph(Surv(survival_time, status) ~ data.matrix(dataframe[, variable]), data=dataframe, na.action=na.omit)
  eta <- fit$linear.predictors
  span <- 1.0 * (length(survival_time[status == 1]) ^ -0.2)
  maintitle <- paste(event, paste(variable, collapse=","),"ROC", sep="_")
  filepath <- paste(basepath, paste(maintitle, extension, sep="."), sep="/")
  OpenOutputdevice(extension, filepath)
  roc <- risksetROC(Stime=survival_time, status=status, marker=eta, predict.time=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  text(0.75, 0.25, paste("AUC =", round(roc$AUC,3)), cex=1.5)
  dev.off()
  maintitle <- paste(event, paste(variable, collapse=","),"AUC", sep="_")
  filepath <- paste(basepath, paste(maintitle, extension, sep="."), sep="/")
  OpenOutputdevice(extension, filepath)
  auc <- risksetAUC(Stime=survival_time, status=status, marker=eta, tmax=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  dev.off()
}

RocAuc <- function(predict_time, survival_time, status, event, variable, dataframe){
  fit <- coxph(Surv(survival_time, status) ~ data.matrix(dataframe[, variable]), data=dataframe, na.action=na.omit)
  eta <- fit$linear.predictors
  span <- 1.0 * (length(survival_time[status == 1]) ^ -0.2)
  maintitle <- paste(event, paste(variable, collapse=","),"ROC", sep="_")
  roc <- risksetROC(Stime=survival_time, status=status, marker=eta, predict.time=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  text(0.75, 0.25, paste("AUC =", round(roc$AUC,3)), cex=1.5)
  dev.copy2pdf(file=paste(basepath, paste(maintitle, extension, sep="."), sep="/"))
  maintitle <- paste(event, paste(variable, collapse=","),"AUC", sep="_")
  auc <- risksetAUC(Stime=survival_time, status=status, marker=eta, tmax=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  dev.copy2pdf(file=paste(basepath, paste(maintitle, extension, sep="."), sep="/"))
}


RocAuc2 <- function(extension, predict_time, survival_time, status, event, variable, dataframe){
  fit <- coxph(Surv(survival_time, status) ~ data.matrix(dataframe[, variable]), data=dataframe, na.action=na.omit)
  eta <- fit$linear.predictors
  span <- 1.0 * (length(survival_time[status == 1]) ^ -0.2)
  maintitle <- paste(event, paste(variable, collapse=","),"ROC", sep="_")
  filepath <- paste(basepath, paste(maintitle, extension, sep="."), sep="/")
  OpenOutputdevice(extension, filepath)
  roc <- risksetROC(Stime=survival_time, status=status, marker=eta, predict.time=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  text(0.75, 0.25, paste("AUC =", round(roc$AUC,3)), cex=1.5)
  dev.off()
  maintitle <- paste(event, paste(variable, collapse=","),"AUC", sep="_")
  filepath <- paste(basepath, paste(maintitle, extension, sep="."), sep="/")
  OpenOutputdevice(extension, filepath)
  auc <- risksetAUC(Stime=survival_time, status=status, marker=eta, tmax=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  dev.off()
}

# Read input data
ads1 <- read.csv('input/170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('input/170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')

# Make ADS
max_therapy_code <- tapply(ads2$REAL_THERAPY_CD, ads2$ID, max)   # extract max therapy code for each patient ID
df1 <- data.frame(ID=row.names(max_therapy_code), max_therapy_code)   # make dataframe only with max therapy code
ads3 <- merge(ads2, df1, by = "ID")
ads6 <- ads3[ads3$REAL_THERAPY_CD == ads3$max_therapy_code & !is.na(ads3$EFF), ]   # only with max therapy code and EFF

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
ads7 <- ads6[ads6$disease_risk == 1 & !is.na(ads6$ARDI2), ]
ads8 <- ads6[ads6$RISK == 1 | ads6$RISK == 2 | ads6$RISK == 3 , ]

basepath <- "/Users/akiko/Dropbox/RDI-Outcome/Figures"
# basepath <- "C:/Users/akiko/Dropbox/RDI-Outcome/Figures"
extension <- "pdf"
RocAuc(1095, ads6$EFS_1, ads6$EFS_FLG, "EFS", c("ARDI", "AGE2C", "disease_risk"), ads6)


# ROC analysis, to see how well the marker predicts 3-year survival
auc <- NULL
main_title <- NULL
predict_time <- NULL
span_ads6_os <- 1.0 * (length(ads6$OS_DAY[ads6$DETH_FLG == 1]) ^ -0.2)
span_ads6_efs <- 1.0 * (length(ads6$EFS_1[ads6$EFS_FLG == 1]) ^ -0.2)
span_ads7_efs <- 1.0 * (length(ads7$EFS_1[ads7$EFS_FLG == 1]) ^ -0.2)
span_ads8_efs <- 1.0 * (length(ads8$EFS_1[ads8$EFS_FLG == 1]) ^ -0.2)
span_ads8_os <- 1.0 * (length(ads8$OS_DAY[ads8$DETH_FLG == 1]) ^ -0.2)

main_title[1] <- "ROC/AUC1: EFS, use ARDI, AGE2c, disease_risk"
predict_time <- 1095
fit1 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta1 <- fit1$linear.predictors
roc1 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta1, predict.time=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads6_efs, main=main_title[1], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc1$AUC,3)), cex=1.5)
auc1 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta1, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_efs, main=main_title[1], lty=2, col="red")

main_title[2] <- "ROC/AUC2: OS, use ARDI, AGE2c, disease_risk"
predict_time <- 1095
fit2 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ ARDI+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta2 <- fit2$linear.predictors
roc2 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta2, predict.time=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[2], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc2$AUC,3)), cex=1.5)
auc2 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta2, tmax=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads6_os, main=main_title[2], lty=2, col="green")

main_title[3] <- "ROC/AUC3: OS, use RDI_ANT, AGE2c, disease_risk"
predict_time <- 1095
fit3 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_ANT+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta3 <- fit3$linear.predictors
roc3 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta3, predict.time=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads6_os, main=main_title[3], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc3$AUC,3)), cex=1.5)
auc3 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta3, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[3], lty=2, col="blue")

main_title[4] <- "ROC/AUC4: OS, use RDI_ARAC, AGE2c, disease_risk"
predict_time <- 1095
fit4 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_ARAC+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta4 <- fit4$linear.predictors
roc4 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta4, predict.time=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads6_os, main=main_title[4], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc4$AUC,3)), cex=1.5)
auc4 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta4, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[4], lty=2, col="pink")

main_title[5] <- "ROC/AUC5: OS, use RDI_VP16, AGE2c, disease_risk"
predict_time <- 1095
fit5 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_VP16+AGE2C+disease_risk, data=ads6, na.action=na.omit)
eta5 <- fit5$linear.predictors
roc5 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta5, predict.time=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads6_os, main=main_title[5], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc5$AUC,3)), cex=1.5)
auc5 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta5, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[5], lty=2, col="purple")

main_title[6] <- "ROC/AUC6: EFS, use ARDI2 and ads7"
predict_time <- 300
fit6 <- coxph(Surv(ads7$EFS_1, ads7$EFS_FLG) ~ ARDI2, data=ads7, na.action=na.omit)
eta6 <- fit6$linear.predictors
roc6 <- risksetROC(Stime=ads7$EFS_1, status=ads7$EFS_FLG, marker=eta6, predict.time=predict_time,
                   method="Schoenfeld",plot=TRUE, span=span_ads7_efs, main=main_title[6], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc6$AUC,3)), cex=1.5)
auc6 <- risksetAUC(Stime=ads7$EFS_1, status=ads7$EFS_FLG, marker=eta6, tmax=predict_time,
                   method="Schoenfeld", span=span_ads7_efs, main=main_title[6], lty=2, col="red")

main_title[7] <- "ROC/AUC7: EFS, use ARDI and ads6"
predict_time <- 1095
fit7 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ ARDI, data=ads6, na.action=na.omit)
eta7 <- fit7$linear.predictors
roc7 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta7, predict.time=predict_time,
                   method="Schoenfeld", plot=TRUE, span=span_ads6_efs, prop=1, main=main_title[7], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc7$AUC,3)), cex=1.5)
auc7 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta7, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_efs, main=main_title[7], lty=2, col="red")

main_title[8] <- "ROC/AUC8: OS, use ARDI1"
predict_time <- 1095
fit8 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ ARDI1, data=ads6, na.action=na.omit)
eta8 <- fit8$linear.predictors
roc8 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta8, predict.time=predict_time,
                   method="Schoenfeld", plot=TRUE, span=span_ads6_os, main=main_title[8], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc8$AUC,3)), cex=1.5)
auc8 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta8, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[8], lty=2, col="red")

main_title[9] <- "ROC/AUC9: EFS, use ARDI1"
predict_time <- 1095
fit9 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ ARDI1, data=ads6, na.action=na.omit)
eta9 <- fit8$linear.predictors
roc9 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta9, predict.time=predict_time,
                   method="Schoenfeld", plot=TRUE, span=span_ads6_efs, main=main_title[9], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc9$AUC,3)), cex=1.5)
auc9 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta9, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_efs, main=main_title[9], lty=2, col="red")

main_title[10] <- "ROC/AUC10: OS, use ARDI2"
predict_time <- 1095
fit10 <- coxph(Surv(ads8$OS_DAY, ads8$DETH_FLG) ~ ARDI2, data=ads8, na.action=na.omit)
eta10 <- fit10$linear.predictors
roc10 <- risksetROC(Stime=ads8$OS_DAY, status=ads8$DETH_FLG, marker=eta10, predict.time=predict_time,
                   method="Schoenfeld", plot=TRUE, span=span_ads8_os, main=main_title[10], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc10$AUC,3)), cex=1.5)
auc10 <- risksetAUC(Stime=ads8$OS_DAY, status=ads8$DETH_FLG, marker=eta10, tmax=predict_time,
                   method="Schoenfeld", span=span_ads8_os, main=main_title[10], lty=2, col="red")

main_title[11] <- "ROC/AUC11: EFS, use ARDI2"
predict_time <- 1095
fit11 <- coxph(Surv(ads8$EFS_1, ads8$EFS_FLG) ~ ARDI2, data=ads8, na.action=na.omit)
eta11 <- fit11$linear.predictors
roc11 <- risksetROC(Stime=ads8$EFS_1, status=ads8$EFS_FLG, marker=eta11, predict.time=predict_time,
                   method="Schoenfeld", plot=TRUE, span=span_ads8_efs, main=main_title[11], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc11$AUC,3)), cex=1.5)
auc11 <- risksetAUC(Stime=ads8$EFS_1, status=ads8$EFS_FLG, marker=eta11, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_efs, main=main_title[11], lty=2, col="red")

main_title[12] <- "ROC/AUC12: OS, use RDI_ARAC RDI_ANT RDI_VP16"
predict_time <- 1095
fit12 <- coxph(Surv(ads6$OS_DAY, ads6$DETH_FLG) ~ RDI_ARAC+RDI_ANT+RDI_VP16, data=ads6, na.action=na.omit)
eta12 <- fit12$linear.predictors
roc12 <- risksetROC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta12, predict.time=predict_time,
                    method="Schoenfeld", plot=TRUE, span=span_ads6_os, main=main_title[12], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc12$AUC,3)), cex=1.5)
auc12 <- risksetAUC(Stime=ads6$OS_DAY, status=ads6$DETH_FLG, marker=eta12, tmax=predict_time,
                   method="Schoenfeld", span=span_ads6_os, main=main_title[12], lty=2, col="red")


main_title[13] <- "ROC/AUC13: EFS, use RDI_ARAC RDI_ANT RDI_VP16"
predict_time <- 1095
fit13 <- coxph(Surv(ads6$EFS_1, ads6$EFS_FLG) ~ RDI_ARAC+RDI_ANT+RDI_VP16, data=ads6, na.action=na.omit)
eta13 <- fit13$linear.predictors
roc13 <- risksetROC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta13, predict.time=predict_time,
                    method="Schoenfeld", plot=TRUE, span=span_ads6_efs, main=main_title[13], lty=2, col="red")
text(0.75, 0.25, paste("AUC =", round(roc12$AUC,3)), cex=1.5)
auc13 <- risksetAUC(Stime=ads6$EFS_1, status=ads6$EFS_FLG, marker=eta13, tmax=predict_time,
                    method="Schoenfeld", span=span_ads6_efs, main=main_title[13], lty=2, col="red")
# cox.zph(fit4)
