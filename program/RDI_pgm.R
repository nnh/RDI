# Read input data
setwd("./input")
ads1 <- read.csv('170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
setwd("..")


survival.time=ads2$EFS_DAY
survival.status=ads2$EFS_FLG

fit0 <- coxph(Surv(EFS_DAY,EFS_FLG)~ ARDI,data=ads2,na.action=na.omit)
eta <- fit0$linear.predictors

ROC.CC10 <- risksetROC(Stime=survival.time,status=survival.status,marker=eta, predict.time=10,method="Cox",data=ads2, na.action=na.omit )
