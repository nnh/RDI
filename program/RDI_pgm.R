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

survival.time=ads4$EFS_DAY
survival.status=ads4$EFS_FLG

fit0 <- coxph(Surv(EFS_DAY,EFS_FLG)~ ARDI,data=ads4,na.action=na.omit)
eta <- fit0$linear.predictors

ROC.CC10 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=10, method="Cox",
                       data=ads4, na.action=na.omit )
