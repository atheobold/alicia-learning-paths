#-----------------------------------------------
#Preliminary data for Fatmeter Calibration
#-----------------------------------------------

install.packages("outliers")

install.packages("lattice")

library(lattice)

#-----------------------------
# WB Lipid Analysis
#-----------------------------

#upper anterior measurement Linear model

str(PADataNoOutlier)

str(PADataNoOutlierMultMeasure)

plot(PADataNoOutlier)

plot(PADataNoOutlierMultMeasure)

linearAnterior <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUA)

summary(linearAnterior)

linearAnterior

with(PADataNoOutlier, plot(Lipid ~ PSUA, las = 1, col = ifelse(PADataNoOutlier$`Fork Length` <  280, "red", "black")))

abline(linearAnterior)

rstudent(linearAnterior)

plot(linearAnterior)

#Exponential function

expAnterior <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUA))

summary (expAnterior)

expAnterior

with(PADataNoOutlier, plot(Lipid ~ log(PSUA), las = 1, 
                           col = ifelse(PADataNoOutlier$`Fork Length` <  260, "red", "black")))

abline(expAnterior)
summary(expAnterior)

plot(expAnterior)

#Cooks Distance for influential point identification

cooksAnterior <- cooks.distance(expAnterior)

plot(cooksAnterior)

abline(h = 4*mean(cooksAnterior, na.rm=T), col="red")

text(x=1:length(cooksAnterior)+1, y=cooksAnterior, labels=ifelse(cooksAnterior>4*mean(cooksAnterior, na.rm=T),names(cooksAnterior),""), col="red")

#Upper posterior measurement

linearposterior <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUP)

summary(linearposterior)

linearposterior

with(PADataNoOutlier, plot(Lipid~ PSUP, las = 1))

abline(linearposterior)

# plot(posterior)

#Exponential posterior measurement 

expPosterior <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUP))

summary(expPosterior)

expPosterior

with(PADataNoOutlier, plot(Lipid ~ log(PSUP)), las = 1)

abline(expPosterior)

# plot(expPosterior)

cooksPosterior <- cooks.distance(expPosterior)

plot(cooksPosterior)

abline(h = 4*mean(cooksPosterior, na.rm=T), col="red")

text(x=1:length(cooksPosterior)+1, y=cooksPosterior, labels=ifelse(cooksPosterior>4*mean(cooksPosterior, na.rm=T),names(cooksPosterior),""), col="red")

grubbs.test(expPosterior)

#outlier test

outlier.test(expPosterior)

#CI

# qt(.975,9)

#upper middle measurements linear

Middle <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM)

summary(Middle)

with(PADataNoOutlier, plot(Lipid ~ PSUM, las = 1))

#upper Middle measurements only 

str(PADataNoOutlier)

expMiddle <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

summary(expMiddle)

with(PADataNoOutlier, plot(Lipid ~ log(PSUM), las = 1, xlab = "Whole-body Lipid Content (%)", 
                           ylab = "UM Fatmeter Reading", col = ifelse(PADataNoOutlier$`Fork Length` <  270, "red", "black")))

abline(expMiddle)

middle

#Cooksd middle

cooksMiddle <- cooks.distance(expMiddle)

plot(cooksMiddle)

abline(h = 4*mean(cooksMiddle, na.rm=T), col="red")

text(x=1:length(cooksMiddle)+1, y=cooksMiddle, labels=ifelse(cooksMiddle>4*mean(cooksMiddle, na.rm=T),names(cooksMiddle),""), col="red")

#---------------------------------
#Energy analysis of data
#---------------------------------

#Means and sd of data

mean(PADataNoOutlier$Energy)

sd(PADataNoOutlier$Energy)

#Anterior energy measurement

expAnteriorE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUA))

summary(expAnteriorE)

with(PADataNoOutlier, plot(Energy ~ log(PSUA), las = 1))

abline(expAnteriorE)

#Posterior energy measurement 

expPosteriorE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUP))

summary(expPosteriorE)

expPosteriorE

with(PADataNoOutlier, plot(Energy ~ log(PSUP), las = 1))

abline(expPosteriorE)

# plot(posteriorE)

#OUTLIER REMOVED anterior Energy

expAnterior2E <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUA))

summary(expAnterior2E)

expAnterior2E

with(PADataNoOutlier, plot(Energy ~ log(PSUA), las = 1))

abline(expAnterior2E)

# plot(anterior2E)

# plot(posterior2E)

# posterior2E

# #CI

# qt(.975,9)

#Middle Data (outlier removed)

expMiddle2E <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM))

summary(expMiddle2E)

expMiddle2E

with(PADataNoOutlier, plot(Energy ~ log(PSUM), las = 1))

abline(expMiddle2E)

# plot(middle2E)

#relationship between lipids and energy

LMLipidEnergy <- lm(PADataNoOutlier$Energy ~ PADataNoOutlier$Lipid)

with(PADataNoOutlier, plot(Energy ~ Lipid, las = 1))

abline(LMLipidEnergy)

summary(LMLipidEnergy)

#relationships with energy and condition

lmConditionEnergy <- lm(PADataNoOutlier$Energy ~ PADataNoOutlier$KN__1)

with(PADataNoOutlier, plot(Energy ~ KN__1, las = 1))

abline(lmConditionEnergy)

summary(lmConditionEnergy)

#relationships with lipids and condition

lmConditionLipid <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$KN__1)

with(PADataNoOutlier, plot(Lipid ~ KN__1, las = 1))

summary(lmConditionLipid)

abline(lmConditionLipid)

#relationship with condition and Fatmeter Anterior

lmConditionAnterior <- lm(PADataNoOutlier$PSUA ~ PADataNoOutlier$KN__1)

with(PADataNoOutlier, plot(PSUA ~ KN__1, las = 1))

abline(lmConditionAnterior)

summary(lmConditionAnterior)

#relationship with condition and Fatmeter Posterior

lmConditionPosterior <- lm(PADataNoOutlier$PSUP ~ PADataNoOutlier$KN__1)

with(PADataNoOutlier, plot(PSUP ~ KN__1, las = 1))

abline (lmConditionPosterior)

summary(lmConditionPosterior)

#relationship with condition and Fatmeter Middle

lmConditionMiddle <- lm(PADataNoOutlier$PSUM ~ PADataNoOutlier$KN__1)

with(PADataNoOutlier, plot(PSUM ~ KN__1, las = 1))

abline(lmConditionMiddle)

summary(lmConditionMiddle)

#combinations

#lm P and A Fatmeter measurements, Kn, and Lipid

lmAPKn <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUA) + log(PADataNoOutlier$PSUP) + PADataNoOutlier$KN__1 + log(PADataNoOutlier$PSUM))

summary(lmAPKn)

plot(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUA) + log(PADataNoOutlier$PSUM) + PADataNoOutlier$KN__1)

abline(lmAPKn)

cookslmAPKn <- cooks.distance(lmAPKn)


plot(cookslmAPKn)

abline(h = 4*mean(cookslmAPKn, na.rm=T), col="red")

text(x=1:length(cookslmAPKn)+1, y=cookslmAPKn, labels=ifelse(cookslmAPKn>4*mean(cookslmAPKn, na.rm=T),names(cookslmAPKn),""), col="red")

#lm of fish middle with length accounted for

lmMidFLM <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM + PADataNoOutlier$'Fork Length' + PADataNoOutlier$Mass +
                 PADataNoOutlier$`Fork Length`*PADataNoOutlier$Mass)
summary(lmMidFLM)

lmMidFL <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM + PADataNoOutlier$'Fork Length')

summary(lmMidFL)

#all size accounted for - Mass not helpful

lmMidSize <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM + PADataNoOutlier$`Fork Length` + 
                  PADataNoOutlier$Mass)

summary(lmMidSize)

#lm of fish middle with mass accounted for

lmMidMass <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM + PADataNoOutlier$Mass)

summary(lmMidMass)

#lm of just fish with multiple measurements

multMiddle <- lm(PADataNoOutlierMultMeasure$Lipid ~ PADataNoOutlierMultMeasure$PSUM)

summary(multMiddle)

plot(PADataNoOutlierMultMeasure$Lipid ~ PADataNoOutlierMultMeasure$PSUM)

abline(multMiddle)

multMiddleExp <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$PSUM))

plot(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$PSUM))

abline(multMiddleExp)

#lm of all measurements

summary(multMiddleExp)

multAll <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$PSUM) + log(PADataNoOutlierMultMeasure$PSUA) + log(PADataNoOutlierMultMeasure$PSUP))

summary(multAll)

plot(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$PSUM) + log(PADataNoOutlierMultMeasure$PSUA) + log(PADataNoOutlierMultMeasure$PSUP))

abline(multAll)

avgAll <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

summary(avgAll)

plot(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

abline(avgAll)

avgAllKn <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$'Avg measurement') + PADataNoOutlierMultMeasure$KN__1)

summary(avgAllKn)

plot(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`) + PADataNoOutlierMultMeasure$KN__1)

abline(avgAllKn)

#alt measurements for energy

lmAPKnE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUA) + log(PADataNoOutlier$PSUP) + PADataNoOutlier$KN__1 + log(PADataNoOutlier$PSUM))

summary(lmAPKnE)

#lm of fish middle with length accounted for

lmMidFLE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length`)

summary(lmMidFLE)

#all size accounted for - Mass not helpful

lmMidSizeE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass)

summary(lmMidSizeE)

#lm of fish middle with mass accounted for

lmMidMassE <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$Mass)

summary(lmMidMassE)

#lm of just fish with multiple measurements

multMiddleE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM))

summary(multMiddleE)

plot(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM))

abline(multMiddleE)

multMiddleExpE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM))

summary(multMiddleExpE)

plot(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM))

abline(multMiddleExpE)

#lm of all measurements

multAllE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM) + log(PADataNoOutlierMultMeasure$PSUA) + log(PADataNoOutlierMultMeasure$PSUP))

summary(multAllE)

plot(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$PSUM) + log(PADataNoOutlierMultMeasure$PSUA) + log(PADataNoOutlierMultMeasure$PSUP))

abline(multAllE)

avgAllE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

summary(avgAllE)

plot(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

abline(avgAllE)

#running Tukey Fences (cutoff 3)

expAnterior

summary(expAnterior)

UMKn <- lm(log)

#LATTICE

library(lattice)

attach(PADataNoOutlier)

Length <- equal.count(PADataNoOutlier$`Fork Length`, number = 3, overlap = 0.1)

str(Length)

densityplot(~Lipid | Length, data = PADataNoOutlier)

bwplot(Lipid~PSUM | Length, data = PADataNoOutlier, layout = c(1,3))

xyplot(Lipid~PSUM | Length, data = PADataNoOutlier, layout = c(1,3))
PADataSmall <- subset(PADataNoOutlier, PADataNoOutlier$'Fork Length' <= 292)

PADataMid <- subset(PADataNoOutlier, PADataNoOutlier$'Fork Length' >= 293 & PADataNoOutlier$'Fork Length' <= 324)

PADataLarge <- subset(PADataNoOutlier, PADataNoOutlier$'Fork Length' >= 325)

Small <- lm(PADataSmall$Lipid ~ PADataSmall$PSUM)

Small

summary(Small)

Mid <- lm(PADataMid$Lipid ~ PADataMid$PSUM)

Mid

summary(Mid)

Large <- lm(PADataLarge$Lipid ~ PADataLarge$PSUM)

Large

summary(Large)

#Favorite Regression

#Middle of all measurements plus Kn

lmKnLipidAllMidKn <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$KN__1)

summary(lmKnLipidAllMidKn) #.2389

#Middle of all measurements plus metric measurements

lmMetricLipidAllMid <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass)
#+PADataNoOutlier$Mass * PADataNoOutlier$`Fork Length`)

summary(lmMetricLipidAllMid) #.2993 (without interaction)//(.4779 with interaction)

plot(lm(PADataNoOutlier$Lipid ~ (log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass
                                 + PADataNoOutlier$Mass * PADataNoOutlier$`Fork Length`))) 

lmMetricLipidInter <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass
                         +PADataNoOutlier$Mass * PADataNoOutlier$`Fork Length`)

lmLipidMid <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

summary(lmLipidMid)

plot(lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))) 

#Middle of all measurements on Lipid

lmAllMid <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

summary(lmAllMid) #.1413

#Avg of P and A

lmPA <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$AvgAP))

summary(lmPA) #.3958

#Avg of AM – good

lmAM <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$AvgAM))


summary(lmAM) #.4606
#Avg of MP

lmMP <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$AvgMP))

summary(lmMP) #.4068

#Avg of AMP - good

lmAMP <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

summary(lmAMP) #.4113

#Condition Factor added to AM and AMP

lmAMKn <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$AvgAM) + PADataNoOutlierMultMeasure$KN)

summary(lmAMKn) #.4545

lmAMPKn <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`) + PADataNoOutlierMultMeasure$KN)

summary(lmAMPKn) #.4257

#Length and Weight added to AM and AMP

lmAMMetric <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$AvgAM) + PADataNoOutlierMultMeasure$`Fork Length`
                 + PADataNoOutlierMultMeasure$Mass + PADataNoOutlierMultMeasure$Mass * PADataNoOutlierMultMeasure$`Fork Length`)

summary(lmAMMetric) #.5318

lmAMPMetric <- lm(PADataNoOutlierMultMeasure$Lipid ~ log(PADataNoOutlierMultMeasure$`Avg measurement`) + PADataNoOutlierMultMeasure$`Fork Length`
                  + PADataNoOutlierMultMeasure$Mass + PADataNoOutlierMultMeasure$Mass * PADataNoOutlierMultMeasure$`Fork Length`)

summary(lmAMPMetric) #.4908

#Avg of P and A

lmPAE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$AvgAP))

summary(lmPAE) #.2714

#Avg of AM - good

lmAME <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$AvgAM))

summary(lmAME) #.3558

#Avg of MP - good

lmMPE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$AvgMP))

summary(lmMPE) #.318

#Avg of AMP

lmAMPE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$`Avg measurement`))

summary(lmAMPE) #.2959

#Condition Factor added to AM and AMP

lmAMKnE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$AvgAM) + PADataNoOutlierMultMeasure$KN)

summary(lmAMKnE) #.3445

lmAMPKnE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$`Avg measurement`) + PADataNoOutlierMultMeasure$KN)

summary(lmAMPKnE) #.3189

#Length and Weight added to AM and AMP

lmAMMetricE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$AvgAM) + PADataNoOutlierMultMeasure$`Fork Length`
                  + PADataNoOutlierMultMeasure$Mass + PADataNoOutlierMultMeasure$Mass * PADataNoOutlierMultMeasure$`Fork Length`)

summary(lmAMMetricE) #.4152

lmAMPMetricE <- lm(PADataNoOutlierMultMeasure$Energy ~ log(PADataNoOutlierMultMeasure$`Avg measurement`) + PADataNoOutlierMultMeasure$`Fork Length`
                   + PADataNoOutlierMultMeasure$Mass + PADataNoOutlierMultMeasure$Mass * PADataNoOutlierMultMeasure$`Fork Length`)

summary(lmAMPMetricE) #.3258

str(PADataNoOutlierMultMeasure)

min(PADataNoOutlierMultMeasure[,3]) #301

min(PADataNoOutlier[,3]) #202

min(PADataNoOutlierMultMeasure[,10], na.rm = T) #.781

min(PADataNoOutlier[,11], na.rm = T) #.796

str(PADataNoOutlier)

#residuals

#Fork  Length

MiddleLM <- lm(Lipid ~ log(PSUM), data = PADataNoOutlier)

str(MiddleLM)

MiddleResid <- resid(MiddleLM)

plot(PADataNoOutlier$`Fork Length`, MiddleResid, ylab = "Residuals", xlab = "Fork Length")

par(mfrow = c(2,2))

#Weight

plot(PADataNoOutlier$Mass, MiddleResid, ylab = "Residuals", xlab = "Mass")

#Kn

plot(PADataNoOutlier$KN__1, MiddleResid, ylab = "Residuals", xlab = "Kn")

#Lipid

plot(PADataNoOutlier$Lipid, MiddleResid, ylab = "Residuals", xlab = "Lipid")


#PSUM

plot(PADataNoOutlier$PSUM, MiddleResid, ylab = "Residuals", xlab = "PSUM")

Lipid <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM)

plot(Lipid)

LipidLog <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

plot(LipidLog)

plot(PADataNoOutlier$Lipid ~ PADataNoOutlier$PSUM)

plot(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

#PSUM + Length

LipidLength <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length`)

summary(LipidLength)

LipidMass <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$Mass)

summary(LipidMass)

anovaLipidAssessment <- anova(lmAllMid, #nocomplex
                              LipidLength, #lowcomplex
                              lmMetricLipidAllMid, #midcomplex
                              lmMetricLipidInter #mostcomplex
)
anovaLipidAssessment

plot(log(PADataNoOutlier$`Fork Length`) ~ log(PADataNoOutlier$Mass))

AIC(LipidLength, LipidMass)

AIC(lmPA, lmAM, lmAMP, lmMP)

anova(lmAM, lmAMMetric)

anova(LipidLength, lmMetricLipidInter)

summary(lmMetricLipidInter)

lmEMKN <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$KN__1)
summary(lmEMKN)

lmEMetricInter <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$Mass 
                     + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass*PADataNoOutlier$`Fork Length`)

summary(lmEMetricInter)

lmEMetric <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$Mass + PADataNoOutlier$`Fork Length`)

summary(lmEMetric)

lmEFL <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length`)

summary(lmEFL)

lmEMass <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$Mass)

summary(lmEMass)

lmEUMAll <- lm(PADataNoOutlier$Energy ~ log(PADataNoOutlier$PSUM))

summary(lmEUMAll)

lmLipidMidInter <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass +PADataNoOutlier$Mass*PADataNoOutlier$`Fork Length`)

summary(lmLipidMidInter)

lmLipidMid <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$`Fork Length`)

summary(lmLipidMid)

anova(lmLipidMid, lmLipidMidInter)

AIC(lmLipidMid, lmLipidMidInter)

Kn <- lm(log(PADataNoOutlier$PSUM) ~ PADataNoOutlier$KN__1)

plot(log(PADataNoOutlier$PSUM) ~ PADataNoOutlier$KN__1)

summary(Kn)

abline(Kn)

lm(log(PADataNoOutlier$Lipid) ~ PADataNoOutlier$KN__1 + PADataNoOutlier$`Fork Length`+ PADataNoOutlier$Mass + PADataNoOutlier$`Fork Length`*PADataNoOutlier$Mass)

summary(lm(log(PADataNoOutlier$Lipid) ~ PADataNoOutlier$KN__1 + PADataNoOutlier$`Fork Length`+ PADataNoOutlier$Mass + PADataNoOutlier$`Fork Length`*PADataNoOutlier$Mass))

lm(log(PADataNoOutlier$Lipid) ~ PADataNoOutlier$PSUM)

summary(lm(log(PADataNoOutlier$Lipid) ~ PADataNoOutlier$PSUM))

plot((PADataNoOutlier$Lipid) ~ log(PADataNoOutlier$`Fork Length`))

plot((PADataNoOutlier$Lipid) ~ log(PADataNoOutlier$Mass))

plot((PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM)))

log <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + log(PADataNoOutlier$Mass) + log(PADataNoOutlier$`Fork Length`) + log(PADataNoOutlier$Mass)*log(PADataNoOutlier$`Fork Length`))

summary(logFL)

summary(log)

logFL <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + log(PADataNoOutlier$`Fork Length`))

anova(log)

anova(lmMetricLipidInter)

summary(lmMetricLipidInter)

anova(lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$`Fork Length` + PADataNoOutlier$Mass))

anova(lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$`Fork Length` * PADataNoOutlier$Mass))

plot(PADataNoOutlier$Lipid ~ PADataNoOutlier$`Fork Length`)

plot(log(PADataNoOutlier$Lipid) ~ log(PADataNoOutlier$`Fork Length`))

KnFat <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM) + PADataNoOutlier$KN__1)

summary(KnFat)

plot(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$PSUM))

plot(KnFat)

Transformed <- lm((PADataNoOutlier$Lipid) ~ log(PADataNoOutlier$PSUM) + log(PADataNoOutlier$`Fork Length`) + log(PADataNoOutlier$Mass) + (log(PADataNoOutlier$`Fork Length`)*log(PADataNoOutlier$Mass)))

summary(Transformed)

plot(Transformed)

LengthMass <- lm((PADataNoOutlier$Lipid) ~ log(PADataNoOutlier$`Fork Length`) + log(PADataNoOutlier$Mass) + (log(PADataNoOutlier$`Fork Length`)*log(PADataNoOutlier$Mass)))

summary(LengthMass)

Kn <- lm(PADataNoOutlier$Lipid ~ PADataNoOutlier$KN__1)

abline(Kn)

plot(PADataNoOutlier$Lipid ~ PADataNoOutlier$Mass)

plot(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$Mass))

plot(PADataNoOutlier$Lipid ~ PADataNoOutlier$`Fork Length`)

plot(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$`Fork Length`))

noInt <- lm(PADataNoOutlier$Lipid ~ log(PADataNoOutlier$`Fork Length`) + log(PADataNoOutlier$Mass))

summary(noInt)

plot(noInt)

large <- PADataNoOutlier[which(PADataNoOutlier$`Fork Length` >= 270 & PADataNoOutlier$'Fork Length' <= 500), ]

small <- PADataNoOutlier[which(PADataNoOutlier$'Fork Length' <= 269), ]

xl <- PADataNoOutlier[which(PADataNoOutlier$`Fork Length`)]

head(large)

plot(large$Lipid ~ large$`Fork Length`)

plot(large$Lipid ~ log(large$`Fork Length`))

largeLM <- lm(large$Lipid ~ large$PSUM + large$`Fork Length` + large$Mass + large$`Fork Length`*large$Mass)

summary(largeLM)

plot(largeLM)

plot(large$Lipid ~ large$Mass)

plot(large$Lipid ~ log(large$Mass))

plot(large$Lipid ~ large$KN__1)

largeLMKn <- lm(large$Lipid ~ large$PSUM + large$KN__1)

summary(largeLMKn)

plot(largeLMKn)

plot(small$Lipid ~ small$`Fork Length`)

plot(small$Lipid ~ small$Mass)

plot(PADataNoPrelim$Lipid ~ PADataNoPrelim$`Fork Length`)

plot(PADataNoPrelim$Lipid ~ log(PADataNoPrelim$`Fork Length`))

plot(PADataNoPrelim$Lipid ~ log(PADataNoPrelim$Mass))

NoPrelimLM <- lm(PADataNoPrelim$Lipid ~ log(PADataNoPrelim$PSUM) + log(PADataNoPrelim$`Fork Length`) + log(PADataNoPrelim$Mass) + log(PADataNoPrelim$`Fork Length`)*log(PADataNoPrelim$Mass))

summary(NoPrelimLM)

plot(PADataNoPrelim$Lipid ~ log(PADataNoPrelim$PSUM))

largeNoPrelim <- PADataNoPrelim[which(PADataNoPrelim$`Fork Length` >= 270 & PADataNoPrelim$'Fork Length' <= 500), ]

smallNoPrelim <- PADataNoPrelim[which(PADataNoPrelim$'Fork Length' <= 269), ]

plot(largeNoPrelim$Lipid ~ log(largeNoPrelim$`Fork Length`))

plot(largeNoPrelim$Lipid ~ log(largeNoPrelim$Mass))

NoPrelimLMPSUM <- lm(PADataNoPrelim$Lipid ~ log(PADataNoPrelim$PSUM))

summary(NoPrelimLMPSUM)

plot(NoPrelimLMPSUM)

NoPrelimBigLM <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`) + log(largeNoPrelim$Mass) + log(largeNoPrelim$`Fork Length`)*log(largeNoPrelim$Mass))

plot(NoPrelimBigLM)

summary(NoPrelimBigLM)

BigPSUMLM <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM))
summary(BigPSUMLM)
BigPUSMLMKn <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$KN...13))
summary(BigPUSMLMKn)  
plot(largeNoPrelim$Lipid ~ log(largeNoPrelim$KN...13))
BigPSUMlmLength <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`))

summary(BigPSUMlmLength)

BigPSUMlmWeight <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$Mass))



summary(BigPSUMlmWeight)
BigPSUMlmLengthWeight <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`) + log(largeNoPrelim$Mass))
summary(BigPSUMlmLengthWeight)

BigPSUMlmLengthWeightInt <- lm(largeNoPrelim$Lipid ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`) + log(largeNoPrelim$Mass) + log(largeNoPrelim$`Fork Length`)*log(largeNoPrelim$Mass))

summary(BigPSUMlmLengthWeightInt)

anovaBigassessment <- anova(BigPSUMLM, #nocomplex
                            BigPSUMlmWeight, #lowcomplex
                            BigPSUMlmLengthWeight, #midcomplex
                            BigPSUMlmLengthWeightInt #mostcomplex
)
summary(anovaBigassessment)

anovaBigassessment

plot(largeNoPrelim$`Fork Length` ~ largeNoPrelim$Mass)

#Energy

BigEnergyLM <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM))

summary(BigEnergyLM)

plot(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM))

BigEnergyLMKn <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$KN...13))

summary(BigPUSMLMKn)  

plot(largeNoPrelim$Energy ~ log(largeNoPrelim$KN...13))

BigEnergylmLength <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`))

summary(BigEnergylmLength)

BigEnergylmWeight <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$Mass))

summary(BigEnergylmWeight)

plot(largeNoPrelim$Energy ~ log(largeNoPrelim$Mass))

BigEnergylmLengthWeight <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`) + log(largeNoPrelim$Mass))

summary(BigEnergylmLengthWeight)

BigEnergylmLengthWeightInt <- lm(largeNoPrelim$Energy ~ log(largeNoPrelim$PSUM) + log(largeNoPrelim$`Fork Length`) + log(largeNoPrelim$Mass) + log(largeNoPrelim$`Fork Length`)*log(largeNoPrelim$Mass))

summary(BigEnergylmLengthWeightInt)

anovaBigEnergyassessment <- anova(BigEnergyLM, #nocomplex
                                  BigEnergylmWeight, #lowcomplex
                                  BigEnergylmLengthWeight, #midcomplex
                                  BigEnergylmLengthWeightInt #mostcomplex
)

summary(anovaBigEnergyassessment)

anovaBigEnergyassessment

plot(largeNoPrelim$`Fork Length` ~ largeNoPrelim $Mass)

install.packages("devtools") 

# required to get packages from GitHub

devtools::install_github("cardiomoon/ggiraphExtra") 

# package from GitHub for interaction plot

library(ggiraphExtra)  

# used to make the interaction plot

largeNoPrelim$
  
  model <- lm(lipid ~ log(fork_length)*log(weight), data = YOUR_DATA) 

model <- lm(Lipid ~ log(`Fork Length`)*log(Mass), data = largeNoPrelim)

# Fit the interaction model you are interested in

ggPredict(model, interactive = TRUE)

# Plots the interaction with different colored points and lines for the weights

ggPredict(BigEnergylmLengthWeightInt, interactive = TRUE)