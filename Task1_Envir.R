# R script: Models in MI Individual Database
# Task 1
# Question 2: effect of environmental variables on sea lice abundance

library(ggplot2)
library(tidyr)
library(corrplot)
library(scales)
library("gridExtra")
library(MASS)
library(dplyr)
library(lme4)
library(AER)
library(DHARMa)
library(vcd)
library(pscl)
library(fitdistrplus)
library(performance)
library(countreg)

# Sea lice database
Copeo<-read.table("Dataset_FIN.csv",
	header = TRUE,sep = ',',fill = TRUE,dec = ".",na.strings = "NA")
str(Copeo)

# Discarding sites with less than 30 samplings and
# fish that were in the seawater more than 30 months
nrow(Copeo[Copeo$TimeSea==30,])
Copeo<-Copeo[Copeo$TimeSea<30,]

Count<-aggregate(Copeo$SampleID,list(Copeo$SiteName),function(x){length(unique(x))})
Discard<-Count[Count$x<30,]$Group.1
Moneo<-Copeo[!(Copeo$SiteName %in% Discard),]

# Discarding covid data (non-integer values)
Moneo<-Moneo[Moneo$Zdroj %in% 'Individ',]

# Making dataset only with variables of interest
Mon<-Moneo[,c(1,3:5,9,19,38,40,56,65:67)]
str(Mon)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of variables
# Question 2: Effect of environmental variables on sea lice abundance 

# On bay level

# Dependent variable: count of total mobile sealice on salmon

# Independent variables:
# bay, temperature, oxygen, salinity, month, year, time in seawater, 
# class of fish, company

# Random effect: nested (Farm/Pen)

# RANDOM EFFECT
# Farm
table(Mon$SiteName)

# Pen
table(Mon$S.R)

# FIXED EFFECT
# Bay
table(Mon$BayName)

# Temperature
summary(Mon$Temp)
hist(Mon$Temp)

# Salinity
summary(Mon$Sal)
hist(Mon$Sal)

# Oxygen
summary(Mon$Oxy1)
hist(Mon$Oxy1)

# Year
unique(Mon$Year)
table(Mon$Year)

# NewMonth
unique(Mon$NewMonth)
table(Mon$NewMonth)

# TimeSea
unique(Mon$TimeSea)
table(Mon$TimeSea)

hist(Moneo$Time)

# Class
unique(Mon$Class)
table(Mon$Class)

# Company
unique(Mon$New_Comp)

Mon[!(Mon$New_Comp %in% c("Mowi","Bifand/Mowi","Bradan Beo Teo.")),]$New_Comp<-'Others'
Mon[Mon$New_Comp %in% c("Bifand/Mowi"),]$New_Comp<-"Mowi"
table(Mon$New_Comp)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Relations between predictors
str(Mon)

# Plots
pairs(Mon[,c(6,2,9:12)],panel=panel.smooth)

# Correlation
res <- cor(Mon[,c(6,2,9:12)])
round(res, 2)

# Changing format
Mon$TimeSea<-as.factor(Mon$TimeSea)
Mon$Year<-as.factor(Mon$Year)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Sea lice DISTRIBUTION

hist(sqrt(Mon$LepeophtheirusSalmonisTotal))

# Proportion of zeros: sea lice total 
T<-nrow(Mon)
Ab<-nrow(Mon[Mon$LepeophtheirusSalmonisTotal==0,])
Ab/T

mean(Mon$LepeophtheirusSalmonisTotal)
var(Mon$LepeophtheirusSalmonisTotal)

# Testing for distribution of sealice total
fit_nb<-goodfit(Mon$LepeophtheirusSalmonisTotal,type='nbinomial') 
fit_po<-goodfit(Mon$LepeophtheirusSalmonisTotal,type='poisson') 
rootogram(fit_nb)
rootogram(fit_po)

# Negative binomial shows better fit than poisson distribution 

distplot(Mon$LepeophtheirusSalmonisTotal, type="poisson")
distplot(Mon$LepeophtheirusSalmonisTotal, type="nbinomial")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# MODELS

# Structure
# LepeophtheirusSalmonisTotal~(1|SiteName/S.R)+BayName+Temp+Sal+Oxy1+Year+NewMonth+
#	TimeSea+Class+New_Comp

# First without random effects

##### Poisson distribution

mod1<-glm(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		NewMonth+TimeSea+Class+New_Comp,
          data=Mon,family='poisson')
summary(mod1)

# Diagnostic plots
par(mfrow=c(2,2))
plot(mod1)

# Testing for overdispersion
E2 <- resid(mod1, type = "pearson")
N  <- nrow(Mon)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)
# Poisson: 10.27

# Checking for zero-inflation
check_zeroinflation(mod1)

##### Quasi-Poisson: correction for overdispersion

mod1.2<-glm(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            family = quasipoisson,
            data = Mon)
summary(mod1.2)

par(mfrow=c(2,2))
plot(mod1.2)

# Checking for zero-inflation
check_zeroinflation(mod1.2)
# Underfitting zeros

#####  Negative binomial distribution

mod2<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
summary(mod2)

# Diagnostic plots
par(mfrow=c(2,2))
plot(mod2)

# Checking for zero-inflation
check_zeroinflation(mod2)
# No reasons for zero-inflation

# Testing for overdispersion

E2 <- resid(mod2, type = "pearson")
N  <- nrow(Mon)
p  <- length(coef(mod2)) + 1  
sum(E2^2) / (N - p)
# Negative binomial: 1.54

# AIC
AIC(mod1,mod2)

# Comparison of Poisson and NB
lrtest(mod1,mod2)
# H0: the Poisson variance equals to the NB variance

# Overall, NB is better fit
# Still problem with overdispersion, but zeros seem to be fine

##### Predictors selection for NB

mod2.1<-glm.nb(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.2<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Sal+Oxy1+Year+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.3<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Oxy1+Year+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.4<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Year+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.5<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+
		NewMonth+TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.6<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		TimeSea+Class+New_Comp,link='log',data=Mon)
mod2.7<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		NewMonth+Class+New_Comp,link='log',data=Mon)
mod2.8<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		NewMonth+TimeSea+New_Comp,link='log',data=Mon)
mod2.9<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
		NewMonth+TimeSea+Class,link='log',data=Mon)

lrtest(mod2,mod2.1);lrtest(mod2,mod2.2)
lrtest(mod2,mod2.3);lrtest(mod2,mod2.4)
lrtest(mod2,mod2.5);lrtest(mod2,mod2.6)
lrtest(mod2,mod2.7);lrtest(mod2,mod2.8)
lrtest(mod2,mod2.9)

# p>0.05= mod2.2, mod2.8
# Discarding: Temp, Class

mod2.10<-glm.nb(LepeophtheirusSalmonisTotal~BayName+Sal+Oxy1+Year+
		NewMonth+TimeSea+New_Comp,link='log',data=Mon)
summary(mod2.10)

# Checking for zero-inflation
check_zeroinflation(mod2.10)

# Testing for overdispersion
E2 <- resid(mod2.10, type = "pearson")
N  <- nrow(Mon)
p  <- length(coef(mod2.10)) + 1  
sum(E2^2) / (N - p)
# Negative binomial: 1.55

AIC(mod2,mod2.1,mod2.2,mod2.3,mod2.5,mod2.6,mod2.7,mod2.8,mod2.9,mod2.10)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##### Zero-altered Poisson (ZAP) vs Zero-altered Negative Binomial (ZANB)

# ZAP
mod3<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'poisson',link='logit',
            data = Mon)
summary(mod3)

# ZANB
mod4<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
summary(mod4)

# Checking for zero-inflation
check_zeroinflation(mod3)
check_zeroinflation(mod4)

# Comparison of the two models
lrtest(mod3,mod4)
AIC(mod3,mod4)

# ZANB is better choice

# Testing for overdispersion
E2 <- resid(mod4, type = "pearson")
N  <- nrow(Mon)
p  <- length(coef(mod4)) + 1
sum(E2^2) / (N - p)
#ZANB: 1.22

# Predictors selection for ZANB
mod4.1<-hurdle(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.2<-hurdle(LepeophtheirusSalmonisTotal~BayName+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.3<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.4<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.5<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.6<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.7<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.8<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
mod4.9<-hurdle(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class,
            dist = 'negbin',link='logit',
            data = Mon)

lrtest(mod4,mod4.1);lrtest(mod4,mod4.2)
lrtest(mod4,mod4.3);lrtest(mod4,mod4.4)
lrtest(mod4,mod4.5);lrtest(mod4,mod4.6)
lrtest(mod4,mod4.7);lrtest(mod4,mod4.8)
lrtest(mod4,mod4.9)

AIC(mod4,mod4.1,mod4.2,mod4.3,mod4.5,mod4.6,mod4.7,mod4.8,mod4.9)

# Model 4, full model is the best

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##### Zero-inflation Poisson (ZIP) vs. Zero-inflation Negative binomial (ZINB)

# All predictors: BayName+S.R+Temp+Sal+Oxy1+Year+NewMonth+TimeSea+Class+New_Com
# Zero-inflation: BayName NewMonth TimeSea 

# ZIP
mod5<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea, 
               dist = 'poisson',
               data = Mon)
summary(mod5)

# ZINB
mod6<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
summary(mod6)

# Comparison of ZIP and ZINB
AIC(mod5,mod6)
lrtest(mod5,mod6)

# Checking for zero-inflation
check_zeroinflation(mod5)
check_zeroinflation(mod6)

# ZINB is better fit

# Testing for overdispersion
E2 <- resid(mod6, type = "pearson")
N  <- nrow(Mon)
p  <- length(coef(mod6))  
sum(E2^2) / (N - p)
# ZINB: 1.57

# Higher overdispersion than with hurdle models

# Predictors selection
mod6.1<-zeroinfl(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.2<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Sal+Oxy1+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.3<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Temp+Oxy1+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.4<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Year+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.5<-zeroinfl(LepeophtheirusSalmonisTotal~BayName+Temp+Sal+Oxy1+NewMonth+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.6<-zeroinfl(LepeophtheirusSalmonisTotal~SiteName+S.R+Temp+Sal+Oxy1+Year+
			TimeSea+Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon) 
mod6.7<-zeroinfl(LepeophtheirusSalmonisTotal~SiteName+S.R+Temp+Sal+Oxy1+Year+NewMonth+
			Class+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.8<-zeroinfl(LepeophtheirusSalmonisTotal~SiteName+S.R+Temp+Sal+Oxy1+Year+NewMonth+
			TimeSea+New_Comp|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)
mod6.9<-zeroinfl(LepeophtheirusSalmonisTotal~SiteName+S.R+Temp+Sal+Oxy1+Year+NewMonth+
			TimeSea+Class|NewMonth+BayName+TimeSea,
                dist='negbin',
                data=Mon)

lrtest(mod6,mod6.1);lrtest(mod6,mod6.2)
lrtest(mod6,mod6.3);lrtest(mod6,mod6.4)
lrtest(mod6,mod6.5);lrtest(mod6,mod6.6)
lrtest(mod6,mod6.7);lrtest(mod6,mod6.9)
# lrtest(mod6,mod6.8)

# Best fit the full model, mod6
AIC(mod6,mod6.1,mod6.2,mod6.3,mod6.4,mod6.5,mod6.6,mod6.7,mod6.8,mod6.9)

# Best fit model without Company, mod6.9

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Comparison of models

# Poisson, NB, ZAP, ZANB, ZIP, ZINB
AIC(mod1,mod2,mod2.10,mod3,mod4,mod5,mod6)

BIC(mod1,mod2,mod2.10,mod3,mod4,mod5,mod6)

Stand_errors<-data.frame(Names=rownames(summary(mod1)$coefficients),
	Poisson=as.numeric(summary(mod1)$coefficients[,2]),
	Quasi.Poisson=as.numeric(summary(mod1.2)$coefficients[,2]),
	NB=as.numeric(summary(mod2)$coefficients[,2]),
	ZAP=as.numeric(summary(mod3)$coefficients$count[,2]),
	ZANB=as.numeric(summary(mod4)$coefficients$count[1:51,2]),
	ZIP=as.numeric(summary(mod5)$coefficients$count[,2]),
	ZINB=as.numeric(summary(mod6)$coefficients$count[1:51,2]))
Stand_errors

apply(Stand_errors[,c(2:8)],2,sum)

# Fit assessed by rootograms
par(mfrow=c(2,3))
rootogram(mod1, max = 200)
rootogram(mod2, max = 200)
rootogram(mod3, max = 200)
rootogram(mod4, max = 200)
rootogram(mod5, max = 200)
rootogram(mod6, max = 200)

# ZANB
rootogram(mod4, main = "ZANB", ylim = c(-5, 15),xlim=c(0,50))
rootogram(mod4, main = "ZANB", ylim = c(-5, 15),xlim=c(20,100))
rootogram(mod4, main = "ZANB", ylim = c(-5, 10),xlim=c(80,200))

rootogram(mod2, main = "NB", ylim = c(-5, 15),xlim=c(0,50))
rootogram(mod2, main = "NB", ylim = c(-5, 15),xlim=c(20,100))
rootogram(mod2, main = "NB", ylim = c(-5, 10),xlim=c(80,200))

# Check Zero-fit of models
check_zeroinflation(mod2)
check_zeroinflation(mod4)
check_zeroinflation(mod6)

# According to AIC the best fit is ZANB, and also the overdispersion
# parameter was the lowest for ZANB
# the standard errors are much higher for all negative binomial distribution
# and by far the highest for ZANB
# Rootograms: best fit is NB, ZANB and ZINB
# Zero-fit: best NB and ZANB

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### Model validation

# Pearson residuals against the fitted values
# To detect non-linearity, unequal error variances and outliers

# YES: residuals randomly distributed around the 0 line
# NO: pattern in the visualization

par(mfrow=c(2,4))
plot(predict(mod1,type='response'),residuals(mod1,type='pearson'),
	main='Poisson', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod1,type='pearson') ~ predict(mod1,type='response')))
plot(predict(mod1.2,type='response'),residuals(mod1.2,type='pearson'),
	main='quasi-Poisson', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod1.2,type='pearson') ~ predict(mod1.2,type='response')))
plot(predict(mod2,type='response'),residuals(mod2,type='pearson'),
	main='NB', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod2,type='pearson') ~ predict(mod2,type='response')))
plot(predict(mod3,type='response'),residuals(mod3,type='pearson'),
	main='ZAP', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod3,type='pearson') ~ predict(mod3,type='response')))
plot(predict(mod4,type='response'),residuals(mod4,type='pearson'),
	main='ZANB', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod4,type='pearson') ~ predict(mod4,type='response')))
plot(predict(mod5,type='response'),residuals(mod5,type='pearson'),
	main='ZIP', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod5,type='pearson') ~ predict(mod5,type='response')))
plot(predict(mod6,type='response'),residuals(mod6,type='pearson'),
	main='ZINB', xlab='Fitted values',ylab='Pearson residuals')
	abline(lm(residuals(mod6,type='pearson') ~ predict(mod6,type='response')))

# The graphs show heteroskelascticity, worse in ~NB models

# Pearson residuals vs. explanatory variables
# ZANB:
# Residuals have no relationship with explanatory variables, 
# they are random values and must not be correlated

EP<-residuals(mod4,type='pearson')
par(mfrow=c(2,5))
plot(Mon$Temp,EP,xlab='Temperature',ylab='Residuals')
	abline(lm(EP ~ Mon$Temp))
plot(Mon$Sal,EP,xlab='Salinity',ylab='Residuals')
	abline(lm(EP ~ Mon$Sal))
plot(Mon$Oxy1,EP,xlab='Oxygen',ylab='Residuals')
	abline(lm(EP ~ Mon$Oxy1))
plot(Mon$Year,EP,xlab='Year',ylab='Residuals')
plot(Mon$TimeSea,EP,xlab='Time in sea',ylab='Residuals')
plot(as.factor(Mon$SiteName),EP,xlab='Site',ylab='Residuals')
plot(as.factor(Mon$NewMonth),EP,xlab='Month',ylab='Residuals')
plot(as.factor(Mon$S.R),EP,xlab='S.R',ylab='Residuals')
plot(as.factor(Mon$Class),EP,xlab='Class',ylab='Residuals')
plot(as.factor(Mon$New_Comp),EP,xlab='Company',ylab='Residuals')

# Should be the same for the whole abline, isnt a big problem

# Observed and fitted values
# OK: all points on a line with an intercept of 0 and a slope of 1 (predicted and observed values
# match up perfectly)
# Points that deviate greatly can indicate outliers or deficiencies in the model

# Poisson
plot(x=predict(mod1), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='Poisson: Predicted vs. Actual Values')
cor(predict(mod1),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3673663
cor(predict(mod1),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5614007

# Quasi-Poisson
plot(x=predict(mod1.2), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='Quasi-Poisson: Predicted vs. Actual Values')
cor(predict(mod1.2),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3669602
cor(predict(mod1.2),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5536822

# NB
plot(x=predict(mod2), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='NB: Predicted vs. Actual Values')
cor(predict(mod2),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3595072
cor(predict(mod2),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.564228

# Prediction of higher values than 4 is not working at all

# ZAP
plot(x=predict(mod3), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='ZAP: Predicted vs. Actual Values')
cor(predict(mod3),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.4830852
cor(predict(mod3),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5726416

# ZANB
plot(x=predict(mod4), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='ZANB: Predicted vs. Actual Values')
cor(predict(mod4),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.4289955
cor(predict(mod4),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5725683

# Prediction is not very accurate, the prediction for high values (more than
# 30 doesnt work at all)

# ZIP
plot(x=predict(mod5), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='ZIP: Predicted vs. Actual Values')
cor(predict(mod5),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 0.4772944
cor(predict(mod5),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5560813

# ZINB
plot(x=predict(mod6), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='ZINB: Predicted vs. Actual Values')
cor(predict(mod6),y=Mon$LepeophtheirusSalmonisTotal,method='pearson')
# 4196446
cor(predict(mod6),y=Mon$LepeophtheirusSalmonisTotal,method='spearman')
# 0.5670201

# Predicted values go up to 40, which is the highest from the models. But observed values reach up to 250

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###### Model assumptions
# For negative binomial hurdle model

### Outliers
plot(x=predict(mod4), y=Mon$LepeophtheirusSalmonisTotal,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='Poisson: Predicted vs. Actual Values')

#### Multicollinearity

check_collinearity(mod4)

# 1=not correlated, 1-5 moderately correlated

# Count: Bayname, Temperature, Year, Month, Company
# Zero: Bayname, Temperature, Month, Company

## Shows high collinearity for several predictors

# Correlation table
Mon$Year<-as.numeric(Mon$Year)
Mon$TimeSea<-as.numeric(Mon$TimeSea)

colnames(Mon)
cor(Mon[,c(2,6,10:13)])

### Temporal Autocorrelation
# Negative Binomial model
# value below two- positive autocorrelation
durbinWatsonTest(residuals(mod2))
# 0.78

ggplot(Mon,aes(as.factor(Year),LepeophtheirusSalmonisTotal))+
	geom_boxplot()
tapply(Mon$LepeophtheirusSalmonisTotal,Mon$Year,mean)

ggplot(Mon,aes(NewMonth,LepeophtheirusSalmonisTotal))+
	geom_boxplot()
tapply(Mon$LepeophtheirusSalmonisTotal,Mon$NewMonth,mean)

md1<-glm.nb(LepeophtheirusSalmonisTotal~Year,link='log',data=Mon)
summary(md1)

par(mfrow=c(2,2))
plot(md1)

# Time-series of residuals
Res<-residuals(md1)
Res1<-sample(Res,50000)

plot(Res1)
plot(Res)
abline(h=0,lty=3)

acf(residuals(md1))

md2<-glm.nb(LepeophtheirusSalmonisTotal~NewMonth,link='log',data=Mon)
summary(md2)

par(mfrow=c(2,2))
plot(md2)

plot(residuals(md2))
abline(h=0,lty=3)

acf(residuals(md1))

# ZANB
durbinWatsonTest(residuals(mod4))
# 1.16: Positive autocorrelation
# Values below 2 mean positive autocorrelation

### Linear relationship between the transformed expected reponse
plot(predict(mod4.2),resid(mod4.2),main='Predicted values vs. Residuals')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##### Focusing on ZANB model

# Random effects as fixed effects
Mon$Year<-as.factor(Mon$Year)
Mon$TimeSea<-as.factor(Mon$TimeSea)

mod4.05<-hurdle(LepeophtheirusSalmonisTotal~BayName+S.R+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)

mod4.06<-hurdle(LepeophtheirusSalmonisTotal~SiteName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)

mod4.07<-hurdle(LepeophtheirusSalmonisTotal~SiteName+S.R+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)

mod4.08<-hurdle(LepeophtheirusSalmonisTotal~S.R+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)

mod4.09<-hurdle(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)

mod4.10<-hurdle(LepeophtheirusSalmonisTotal~BayName+SiteName+Temp+Sal+Oxy1+Year+
			NewMonth+TimeSea+Class+New_Comp,
            dist = 'negbin',link='logit',
            data = Mon)
summary(mod4.10) # some zero hurdle model coefficients are aliased, regressors dropped 


AIC(mod4,mod4.05,mod4.06,mod4.07,mod4.08,mod4.09,mod4.10)
AIC(mod4,mod4.05,mod4.06,mod4.07,mod4.08,mod4.09)

AIC(mod4.06,mod4.07)
### mod4.06 and mod4.07 the best AIC values
### Site instead of Bay, with or without Pen

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




