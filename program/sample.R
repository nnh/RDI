library(survival)
library(risksetROC)
library(MASS)
data(VA)
survival.time <- VA$stime
survival.status <- VA$status
score <- VA$Karn
cell.type <- factor(VA$cell)
tx <- as.integer(VA$treat == 1)
age <- VA$age
survival.status[survival.time > 500 ] <- 0
survival.time[survival.time > 500 ] <- 500
fit0 <- coxph(Surv(survival.time, survival.status) ~ score + cell.type + tx + age, na.action=na.omit)
summary(fit0)
eta <- fit0$linear.predictor
ROC.CC30 <- risksetROC(Stime=survival.time, status=survival.status, marker=eta, predict.time=30,
                       method="Cox", main="ROC Curve", lty=2, col="red")
AUC <- NULL
out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status, predict.time=30)
AUC <- out$AUC   ## to see how well the marker predicts one-month survival
