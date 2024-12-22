names(RPMA2Growth)

early <- subset(RPMA2Growth, StockYear < 2006)

mid <- subset(RPMA2Growth, StockYear < 2014 & StockYear > 2003)

RPMA2GrowthSub <- transform(RPMA2Growth, Age = as.integer(Age))

Early <- subset(RPMA2GrowthSub, StockYear < 2004)

Mid <- subset(RPMA2GrowthSub, StockYear < 2018 & StockYear > 2005)

EarlyWeightAge <- ddply(Early, ~Age, summarise, meanWE=mean(Weight, na.rm = T))

EarlyLengthAge <- ddply(Early, ~Age, summarise, meanLE = mean(ForkLength, na.rm = T))

MidLengthAge <- ddply(Mid, ~Age, summarise, meanLM = mean(ForkLength, na.rm = T))

WeightChange <- rep(NA, 9)

library(plyr)

WeightAge <- ddply(RPMA2GrowthSub, ~Age, summarise, meanW=mean(Weight, na.rm = T))

LengthAge <- ddply(RPMA2GrowthSub, ~Age, summarise, meanL=mean(ForkLength, na.rm = T))

#Tanner's code/help

WeightChange <- rep(NA, 9)

library(plyr)

WeightAge <- ddply(RPMA2GrowthSub, ~Age, summarise, meanW=mean(Weight, na.rm = T))

LengthAge <- ddply(RPMA2GrowthSub, ~Age, summarise, meanL=mean(ForkLength, na.rm = T))

plot(WeightAge$meanW ~ WeightAge$Age)

plot(LengthAge$mean ~ LengthAge$Age)

WeightChange

Weight1 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 1], na.rm = TRUE)

Weight1

Length1 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 1], na.rm  = TRUE)

Weight2 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 2], na.rm = TRUE)

Length2 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 2], na.rm  = TRUE)

Weight3 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 3], na.rm = TRUE)

Length3 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 3], na.rm  = TRUE)

Weight4 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 4], na.rm = TRUE)

Length4 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 4], na.rm  = TRUE)

Weight5 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 5], na.rm = TRUE)

Length5 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 5], na.rm  = TRUE)

Weight6 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 6], na.rm = TRUE)

Length6 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 6], na.rm  = TRUE)

Weight7 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 7], na.rm = TRUE)

Length7 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 7], na.rm  = TRUE)

Weight8 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 8], na.rm = TRUE)

Length8 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 8], na.rm  = TRUE)

Weight9 <- mean(RPMA2GrowthSub$Weight[RPMA2GrowthSub$Age == 9], na.rm = TRUE)

Length9 <- mean(RPMA2GrowthSub$ForkLength[RPMA2GrowthSub$Age == 9], na.rm  = TRUE)

