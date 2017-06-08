# Read input data
setwd("./input")
ads1 <- read.csv('170427_ADS1.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
ads2 <- read.csv('170427_ADS2.csv', as.is = T, fileEncoding = 'UTF-8-BOM')
setwd("..")