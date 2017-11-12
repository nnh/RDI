library(survival)
library(risksetROC)

RocAuc <- function(predict_time, survival_time, status, event, variable, dataframe){
  fit <- coxph(Surv(survival_time, status) ~ data.matrix(dataframe[, variable]), data=dataframe, na.action=na.omit)
  eta <- fit$linear.predictors
  span <- 1.0 * (length(survival_time[status == 1]) ^ -0.2)
  maintitle <- paste(event, paste(variable, collapse=","),"ROC", sep="_")
  roc <- risksetROC(Stime=survival_time, status=status, marker=eta, predict.time=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  text(0.75, 0.25, paste("AUC =", round(roc$AUC, 3)), cex=1.5)
  # dev.copy2pdf(file=paste(basepath, paste(maintitle, "pdf", sep="."), sep="/"))
  dev.copy(png, file=paste(basepath, paste(maintitle, "png", sep="."), sep="/"), width=1000)
  dev.off()
  maintitle <- paste(event, paste(variable, collapse=","),"AUC", sep="_")
  auc <- risksetAUC(Stime=survival_time, status=status, marker=eta, tmax=predict_time,
                    method="Schoenfeld", span=span, main=maintitle, lty=2, col="red")
  # dev.copy2pdf(file=paste(basepath, paste(maintitle, "pdf", sep="."), sep="/"))
  dev.copy(png, file=paste(basepath, paste(maintitle, "png", sep="."), sep="/"), width=1000)
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

# basepath <- "/Users/akiko/Dropbox/RDI-Outcome/Figures"
basepath <- "C:/Users/akiko/Dropbox/RDI-Outcome/Figures"

RocAuc(1095, ads6$EFS_1, ads6$EFS_FLG, "EFS", c("ARDI", "AGE2C", "disease_risk"), ads6)
RocAuc(1095, ads6$EFS_1, ads6$EFS_FLG, "EFS", c("AGE2C", "disease_risk"), ads6)
RocAuc(1095, ads6$EFS_1, ads6$EFS_FLG, "EFS", c("ARDI"), ads6)
RocAuc(1095, ads6$EFS_1, ads6$EFS_FLG, "EFS", c("ARDI1"), ads6)
RocAuc(1095, ads8$EFS_1, ads8$EFS_FLG, "EFS", c("ARDI2"), ads8)

RocAuc(1095, ads6$OS_DAY, ads6$DETH_FLG, "OS", c("ARDI", "AGE2C", "disease_risk"), ads6)
RocAuc(1095, ads6$OS_DAY, ads6$DETH_FLG, "OS", c("AGE2C", "disease_risk"), ads6)
RocAuc(1095, ads6$OS_DAY, ads6$DETH_FLG, "OS", c("ARDI"), ads6)
RocAuc(1095, ads6$OS_DAY, ads6$DETH_FLG, "OS", c("ARDI1"), ads6)
RocAuc(1095, ads8$OS_DAY, ads8$DETH_FLG, "OS", c("ARDI2"), ads8)
