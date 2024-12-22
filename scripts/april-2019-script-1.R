#Preliminary data for Fatmeter Calibration
#-----------------------------------------------

#-----------------------------
# WB Lipid Analysis
#-----------------------------

#upper anterior measurement(outlier included) Linear model

#Only fish that had upper anterior measurements were the larger cohort sampled at Garrison

linearAnterior <- lm(ProximateAnalysisData$PSUA~ProximateAnalysisData$Lipid)

summary(linearAnterior)

linearAnterior

with(ProximateAnalysisData, plot(PSUA ~ Lipid, las = 1))

abline(linearAnterior)

plot(linearAnterior)

#Exponential function

expAnterior <- lm(log(ProximateAnalysisData$PSUA)~ProximateAnalysisData$Lipid)

summary(expAnterior)

expAnterior

with(ProximateAnalysisData, plot(log(PSUA) ~ Lipid, las = 1))

#Upper posterior measurement (outlier included)

#Only fish that had upper posterior measurements were the larger cohort sampled at Garrison

linearposterior <- lm(ProximateAnalysisData$PSUP ~ ProximateAnalysisData$Lipid)

summary(linearposterior)

linearposterior
with(ProximateAnalysisData, plot(PSUP ~ Lipid, las = 1))

plot(posterior)

#Exponential posterior measurement (outlier included)

expPosterior <- lm(log(ProximateAnalysisData$PSUP) ~ ProximateAnalysisData$Lipid)

summary(expPosterior)

expPosterior

with(ProximateAnalysisData, plot(log(PSUP) ~ Lipid), las = 1)

plot(expPosterior)

#Anterior measuremnet with OUTLIER REMOVED exponential

expAnterior2 <- lm(log(ProximateAnalysisDataOutlier$PSUA) ~ ProximateAnalysisDataOutlier$Lipid)

summary(expAnterior2)

expAnterior2

with(ProximateAnalysisDataOutlier, plot(log(PSUA) ~ Lipid), las = 1,  xlab = "Whole-body Lipid Content (%)", ylab = "UA Fatmeter Reading")

abline(expAnterior2)

plot(anterior2)

anterior2

#Posterior measuremnet with OUTLIER REMOVED

expPosterior2 <- lm(log(ProximateAnalysisDataOutlier$PSUP) ~ ProximateAnalysisDataOutlier$Lipid)

summary(expPosterior2)

with(ProximateAnalysisDataOutlier, plot(log(PSUP) ~ Lipid, las = 1))

abline(expPosterior2)

plot(posterior2)

posterior2

#CI

qt(.975,9)

#upper Middle measurements only (Including outlier)

expMiddle <- lm(ProximateAnalysisData$PSUM ~ log(ProximateAnalysisData$Lipid))

summary(expMiddle)

with(ProximateAnalysisData, plot(PSUM ~ log(Lipid), las = 1, xlab = "Whole-body Lipid Content (%)", ylab = "UM Fatmeter Reading"))

plot(middle)

middle

#---------------------------------
#Energy analysis of data
#---------------------------------

#Means and sd of data

mean(ProximateAnalysisData$Energy)

sd(ProximateAnalysisData$Energy)

#Anterior energy measurement (outlier included)

expAnteriorE <- lm(log(ProximateAnalysisData$PSUA)~ProximateAnalysisData$Energy)

summary(expAnteriorE)

with(ProximateAnalysisData, plot(log(PSUA) ~ Energy, las = 1))

plot(anteriorE)

#Posterior energy measurement (outlier included)

expPosteriorE <- lm(log(ProximateAnalysisData$PSUP) ~ ProximateAnalysisData$Energy)

summary(expPosteriorE)

expPosteriorE

with (ProximateAnalysisData, plot(PSUP ~ Energy, las = 1))

plot(posteriorE)

#OUTLIER REMOVED anterior Energy

expAnterior2E <- lm(log(ProximateAnalysisDataOutlier$PSUA) ~ ProximateAnalysisDataOutlier$Energy)

summary(expAnterior2E)

expAnterior2E

with(ProximateAnalysisDataOutlier, plot(log(PSUA) ~ Energy, las = 1))

abline(expAnterior2E)

plot(anterior2E)

#Posterior energy (outlier removed)

expPosterior2E <- lm(log(ProximateAnalysisDataOutlier$PSUP) ~ ProximateAnalysisDataOutlier$Energy)

summary(expPosterior2E)

with(ProximateAnalysisDataOutlier, plot(log(PSUP) ~ Energy, las = 1))

abline(expPosterior2E)

plot(posterior2E)

posterior2E

#CI

qt(.975,9)

#Middle Data (outlier removed)

expMiddle2E <- lm(log(ProximateAnalysisDataOutlier$PSUM) ~ ProximateAnalysisDataOutlier$Energy)

summary(expMiddle2E)

expMiddle2E

with(ProximateAnalysisDataOutlier, plot(PSUM ~ log(Energy), las = 1))

plot(middle2E)