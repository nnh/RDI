library(survival)
library(risksetROC)
# Read input data
ads1 <- read.csv('input/170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('input/170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')

MAX_T <- tapply(ads2$REAL_THERAPY_CD, ads2$ID, max)
max_therapy <- data.frame(MAX_T)
max_therapy$ID <- row.names(MAX_T)
ads3 <- merge(ads2, max_therapy, by = "ID")
ads4 <- ads3[ads3$REAL_THERAPY_CD == ads3$MAX_T, ]
ads5 <- ads4[!is.na(ads4$EFS_DAY) & !is.na(ads4$EFS_FLG), ]

survival.time <- ads5$EFS_DAY
survival.status <- ads5$EFS_FLG

fit0 <- coxph(Surv(EFS_DAY, EFS_FLG) ~ ARDI, data=ads5, na.action=na.omit)
eta <- fit0$linear.predictors
ROC.CC10 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=10, method="Cox",
                       main="ROC Curve", lty=2, col="red")
AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status, predict.time=10)
AUC <- out$AUC

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

