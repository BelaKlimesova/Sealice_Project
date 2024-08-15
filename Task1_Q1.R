
# Task 1
# Question 1: The effect of environmental variables on sea lice abundance

library(ggplot2)
library(glmmTMB)
library(xts)
library(car)
library(lmtest)
library(gridExtra)

# Sea lice dataset
Copeo<-read.table("Dataset_FIN.csv",
	header = TRUE,sep = ',',fill = TRUE,dec = ".",na.strings = "NA")
str(Copeo)

# Discarding sites with less than 30 samplings and
# fish that were in the seawater more than 30 months
nrow(Copeo[Copeo$TimeSea==30,])
Copeo<-Copeo[Copeo$TimeSea<30,]

Pocty<-aggregate(Copeo$SampleID,list(Copeo$SiteName),function(x){length(unique(x))})
Jmen_vyr<-Pocty[Pocty$x<30,]$Group.1
Moneo<-Copeo[!(Copeo$SiteName %in% Jmen_vyr),]

# Discarding covid data (non-integer values)
Mon<-Moneo[Moneo$Zdroj %in% 'Individ',]

# Changing format of variables
Mon$Year_f<-as.factor(Mon$Year)
Mon$TimeSea_f<-as.factor(Mon$TimeSea)

Mon[!(Mon$New_Comp %in% c("Mowi","Bifand/Mowi","Bradan Beo Teo.")),]$New_Comp<-'Others'
Mon[Mon$New_Comp %in% c("Bifand/Mowi"),]$New_Comp<-"Mowi"

# Discarding outliers
Mon_2<-Mon[Mon$LepeophtheirusSalmonisTotal<161,]
range(Mon_2$LepeophtheirusSalmonisTotal)

Mon_2[c(165496), ]
Mon_2.1<- Mon_2[-c(165496), ]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################################################################################################################

############################# Fish level

# New variables
Mon_2.1$Year_n<-as.numeric(Mon_2.1$Year_f)/10
Mon_2.1$Year_w<-numFactor(Mon_2.1$Year)
Mon_2.1$Month_w<-numFactor(Mon_2.1$Month)
Mon_2.1$Site<-as.numeric(as.factor(Mon_2.1$SiteName))

# Exploration of data autocorrelation ACF
Mon_2.1$Date<-as.Date(Mon_2.1$Date)   
Mmm<-xts(Mon_2.1$LepeophtheirusSalmonisTotal,Mon_2.1$Date) 

acf(coredata(Mmm),1000)

# Model without covariate structure for temporal autocorrelation
m1<-glmmTMB(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year_n+
		TimeSea_f+Class+New_Comp+S.R+(1|SiteName),
		family=nbinom2,data=Mon_2.1)
summary(m1)

# Model with covariate structure: Year_w
m1.1<-glmmTMB(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year_n+
		TimeSea_f+Class+New_Comp+S.R+(1|SiteName)+ou(Year_w+0|Site),
		family=nbinom2,data=Mon_2.1)
summary(m1.1)

# Model with covariate structure: Month_w
m1.2<-glmmTMB(LepeophtheirusSalmonisTotal~Temp+Sal+Oxy1+Year_n+
		TimeSea_f+Class+New_Comp+S.R+(1|SiteName)+ou(Month_w+0|Site),
		family=nbinom2,data=Mon_2.1)
summary(m1.2)

# Checking autocorrelation and comparison of three models
durbinWatsonTest(resid(m1))
durbinWatsonTest(resid(m1.1))
durbinWatsonTest(resid(m1.2))

dwtest(m1)
dwtest(m1.1)
dwtest(m1.2)

bgtest(m1)
bgtest(m1.1)
bgtest(m1.2)

acf(resid(m1),1000)
acf(resid(m1.1),1000)
acf(resid(m1.2),1000)

# AIC
AIC(m1,m1.1,m1.2)

# Correlation of predicted values and raw data
cor(predict(m1,type="response"),y=Mon_2.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3897307
cor(predict(m1.1,type="response"),y=Mon_2.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.5691623
cor(predict(m1.2,type="response"),y=Mon_2.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.4938102

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results

# only looking at variables of interest: temperature, salinity, oxygen
# using m1: without covariate structure

# Temperature

newdata2<- data.frame(
	Temp = rep(seq(from = min(Mon_2.1$Temp), to = max(Mon_2.1$Temp),length.out = 100),4),
	Sal = mean(Mon_2.1$Sal),Oxy1= mean(Mon_2.1$Oxy1),Year_n=mean(Mon_2.1$Year_n),
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_2.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

A<-ggplot(newdata2, aes(Temp,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Temperature", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Salinity
newdata2<- data.frame(
	Temp = mean(Mon_2.1$Temp),
	Sal = rep(seq(from = min(Mon_2.1$Sal), to = max(Mon_2.1$Sal),length.out = 100),4),
	Oxy1= mean(Mon_2.1$Oxy1),Year_n=mean(Mon_2.1$Year_n),
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_2.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

B<-ggplot(newdata2, aes(Sal,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Salinity", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Oxygen
newdata2<- data.frame(
	Temp = mean(Mon_2.1$Temp),
	Sal = mean(Mon_2.1$Sal),
	Oxy1= rep(seq(from = min(Mon_2.1$Oxy1), to = max(Mon_2.1$Oxy1),length.out = 100),4),
	Year_n=mean(Mon_2.1$Year_n),
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_2.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Others',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

C<-ggplot(newdata2, aes(Oxy1,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Oxygen", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(A,B,C,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

############################################################################################################################################

############################# Pen level

# New dataset on the level of pen: mean abundance per pen

Mon<-aggregate(Mon_2.1$LepeophtheirusSalmonisTotal,list(
	Mon_2.1$Year_n,Mon_2.1$Year,Mon_2.1$SiteName,Mon_2.1$Month,
	Mon_2.1$NewMonth,Mon_2.1$S.R,Mon_2.1$Class,Mon_2.1$TimeSea,Mon_2.1$Temp,
	Mon_2.1$Sal,Mon_2.1$Oxy1,Mon_2.1$New_Comp,Mon_2.1$SampleID),mean)
colnames(Mon)<-c('Year','Year_f','SiteName','Month','NewMonth','S.R','Class','TimeSea',
	'Temp','Sal','Oxy','Company','ID','Lep')
head(Mon)

# New variables
Mon$TimeSea<-as.factor(Mon$TimeSea)

# Exploration of autocorrelation for data
Mon$Date<-paste0(Mon$Year_f,'/',Mon$Month,'/01')
Mon$Date<-as.Date(Mon$Date)
Mmm<-xts(Mon$Lep,Mon$Date) 

acf(coredata(Mmm),1000)

# Full Model
m2<-glmmTMB(Lep~Temp+Sal+Oxy+Year+TimeSea+Class+
		Company+S.R+(1|SiteName),
		ziformula=~Temp+Sal+Oxy+Year+TimeSea+Class+Company+S.R,
		family=ziGamma(link="log"),data=Mon)
summary(m2)

# Checking autocorrelation
durbinWatsonTest(resid(m2))
dwtest(m2)
bgtest(m2)

acf(resid(m2),1000)

# Correlation of predicted values and reponse variable
cor(predict(m2,type="response"),y=Mon$Lep,method='pearson')
# 0.4383131

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results for pen level
# Only variables of interest
summary(m2)

# Temperature

newdata2<- data.frame(
	Temp = rep(seq(from = min(Mon$Temp), to = max(Mon$Temp),length.out = 100),4),
	Sal = mean(Mon$Sal),Oxy= mean(Mon$Oxy),Year=mean(Mon$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

D<-ggplot(newdata2, aes(Temp,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Temperature", y = "Abundance")+ 
	ylab('')+
  ggtitle('Pen level')+
  theme(legend.position = "none",text = element_text(size = 17))

# Salinity
newdata2<- data.frame(
	Temp = mean(Mon$Temp),
	Sal = rep(seq(from = min(Mon$Sal), to = max(Mon$Sal),length.out = 100),4),
	Oxy= mean(Mon$Oxy),Year=mean(Mon$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

E<-ggplot(newdata2, aes(Sal,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Salinity", y = "Abundance")+ 
	ylab('')+
  ggtitle('Pen level')+
  theme(legend.position = "none",text = element_text(size = 17))

# Oxygen
newdata2<- data.frame(
	Temp = mean(Mon$Temp),Sal = mean(Mon$Sal),
	Oxy= rep(seq(from = min(Mon$Oxy), to = max(Mon$Oxy),length.out = 100),4),
	Year=mean(Mon$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),S.R =rep('Random',400)
	)

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

F<-ggplot(newdata2, aes(Oxy,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Oxygen", y = "Abundance")+ 
	ylab('')+
  ggtitle('Pen level')+
  theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(D,E,F,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################################################################################################################

############################# Farm level

# New dataset on the level of farm: mean per farm

Moi<-aggregate(Mon_2.1$LepeophtheirusSalmonisTotal,list(
	Mon_2.1$Year_n,Mon_2.1$Year,Mon_2.1$SiteName,Mon_2.1$NewMonth,Mon_2.1$Month,
	Mon_2.1$Class,Mon_2.1$TimeSea,Mon_2.1$Date,Mon_2.1$Temp,
	Mon_2.1$Sal,Mon_2.1$Oxy1,Mon_2.1$New_Comp),mean)
colnames(Moi)<-c('Year','Year_f','SiteName','NewMonth','Month','Class',
	'TimeSea','Date','Temp','Sal','Oxy','Company','Lep')
head(Moi)

# New variables
Moi$TimeSea<-as.factor(Moi$TimeSea)

# Exploration of autocorrelation for data
Moi$Date<-paste0(Moi$Year_f,'/',Moi$Month,'/01')
Moi$Date<-as.Date(Moi$Date)
Mmm<-xts(Moi$Lep,Moi$Date) 

acf(coredata(Mmm),1000)

# Full Model
m3<-glmmTMB(Lep~Temp+Sal+Oxy+Year+TimeSea+Class+Company+(1|SiteName),
		ziformula=~Temp+Sal+Oxy+Year+TimeSea+Class+Company,
		family=ziGamma(link="log"),data=Moi)
summary(m3)

# Checking autocorrelation and comparison of models
durbinWatsonTest(resid(m3))
dwtest(m3)
bgtest(m3)

acf(resid(m3),1000)

# Correlation of predicted values and reponse variable
cor(predict(m3,type="response"),y=Moi$Lep,method='pearson')
# 0.455614

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results for farm level
# only variables of interest
summary(m3)

# Temperature
newdata2<- data.frame(
	Temp = rep(seq(from = min(Moi$Temp), to = max(Moi$Temp),length.out = 100),4),
	Sal = mean(Moi$Sal),Oxy= mean(Moi$Oxy),Year=mean(Moi$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Moi$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

G<-ggplot(newdata2, aes(Temp,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Temperature", y = "Abundance")+ 
	ylab('')+
  ggtitle('Farm level')+
  theme(legend.position = "none",text = element_text(size = 17))

# Salinity
newdata2<- data.frame(
	Temp = mean(Moi$Temp),
	Sal = rep(seq(from = min(Moi$Sal), to = max(Moi$Sal),length.out = 100),4),
	Oxy= mean(Moi$Oxy),Year=mean(Moi$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Moi$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

H<-ggplot(newdata2, aes(Sal,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Salinity", y = "Abundance")+ 
	ylab('')+
  ggtitle('Farm level')+
  theme(legend.position = "none",text = element_text(size = 17))

# Oxygen
newdata2<- data.frame(
	Temp = mean(Moi$Temp),Sal = mean(Moi$Sal),
	Oxy= rep(seq(from = min(Moi$Oxy), to = max(Moi$Oxy),length.out = 100),4),
	Year=mean(Moi$Year),
	TimeSea = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Moi$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

I<-ggplot(newdata2, aes(Oxy,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Oxygen", y = "Abundance")+ 
	ylab('')+
  ggtitle('Farm level')+
  theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(G,H,I,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################################################################################################

# Comparison of variation between datasets

# Fish level
sd(Mon_2.1$LepeophtheirusSalmonisTotal)
mean(Mon_2.1$LepeophtheirusSalmonisTotal)
median(Mon_2.1$LepeophtheirusSalmonisTotal)
summary(Mon_2.1$LepeophtheirusSalmonisTotal)

# Pen level
sd(Mon$Lep)
mean(Mon$Lep)
median(Mon$Lep)
summary(Mon$Lep)

# Farm level
sd(Moi$Lep)
mean(Moi$Lep)
median(Moi$Lep)
summary(Moi$Lep)

# Number of observation per dataset
# Fish level
nrow(Mon_2.1)
# Pen level
nrow(Mon)
# Farm level
nrow(Moi)

# Comparison of models

# Fish level
summary(m1)

# Pen level
summary(m2)

# Farm level
summary(m3)

# Comparison of results
grid.arrange(A,D,G,ncol = 3)
grid.arrange(B,E,H,ncol = 3)
grid.arrange(C,F,I,ncol = 3)

# Everything together
grid.arrange(A,B,C,D,E,F,G,H,I,ncol = 3)



