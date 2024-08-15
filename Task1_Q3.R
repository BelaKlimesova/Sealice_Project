
# Task 1: Epidemiology study
# Question 3: Long term changes in sea lice abundance in Ireland (1995-2022) 

library(ggplot2)
library(glmmTMB)
library(xts)
library(car)
library(lmtest)
library(gridExtra)

# Sea lice database
Meany<-read.table("Mean_FIN.csv",
	header = TRUE,sep = ',',fill = TRUE,dec = ".",na.strings = "NA")

# Discarding sites and fish long at sea
Pocty<-aggregate(Meany$SplDate,list(Meany$SiteName),function(x){length(unique(x))})
Jmen_1rok<-Pocty[Pocty$x<120,]$Group.1

Coli<-Meany[Meany$TimeSea<30,]
Copeo<-Coli[!(Coli$SiteName %in% Jmen_1rok),]

Copeo$SplDate<-as.Date(Copeo$SplDate)

# Making dataset only with variables of interest
# Discarding data higher than 90
Mon<-Copeo[Copeo$MeanLepTotal<91,]

# Adjustment of variables
Mon$Year_f<-as.factor(Mon$Year)
Mon$TimeSea_f<-as.factor(Mon$TimeSea)
Mon$CageType<-as.factor(Mon$CageType)
Mon$BayName<-as.factor(Mon$BayName)
Mon$Class<-as.factor(Mon$Class)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################################################################################################################

############################# Pen level

# New variables
Mon$Year_n<-as.numeric(Mon$Year_f)/10
Mon$Year_w<-numFactor(Mon$Year)
Mon$Site<-as.numeric(as.factor(Mon$SiteName))

Mon$NewMonth<-as.factor(Mon$NewMonth)
Mon$NewMonth<-factor(Mon$NewMonth,levels=c('Jan/Dec','Feb','Mar1','Mar2',
	'Apr1','Apr2','May1','May2','Jun','Jul','Aug','Sep','Oct','Nov'))
Mon$Month_w<-numFactor(as.numeric(Mon$NewMonth))

# Exploration of autocorrelation of data
Mon$Date<-as.Date(Mon$SplDate)   
Mmm<-xts(Mon$MeanLepTotal,Mon$Date) 

acf(coredata(Mmm),1000)

# Model without a covariate structure for temporal autocorrelation
m1.0<-glmmTMB(MeanLepTotal~Class+Year_f+NewMonth+
	TimeSea_f+CageType+BayName+(1|SiteName),
	ziformula=~Class+Year_f+NewMonth+CageType+BayName+TimeSea_f,
	family=ziGamma(link="log"),data=Mon)
summary(m1.0)

# Model with covariate structure: Year_w
m1.1<-glmmTMB(MeanLepTotal~Class+Year_f+NewMonth+
	TimeSea_f+CageType+BayName+(1|SiteName)+ou(Year_w+0|Site),
	ziformula=~Class+Year_f+NewMonth+CageType+BayName+TimeSea_f,
	family=ziGamma(link="log"),data=Mon)
summary(m1.1)

# Model with covariate structure: Month_w
m1.2<-glmmTMB(MeanLepTotal~Class+Year_f+NewMonth+
	TimeSea_f+CageType+BayName+(1|SiteName)+ou(Month_w+0|Site),
	ziformula=~Class+Year_f+NewMonth+CageType+BayName+TimeSea_f,
	family=ziGamma(link="log"),data=Mon)
summary(m1.2)

# Checking autocorrelation and comparison of three models
durbinWatsonTest(resid(m1.0))
durbinWatsonTest(resid(m1.1))
durbinWatsonTest(resid(m1.2))

dwtest(m1.0)
dwtest(m1.1)
dwtest(m1.2)

bgtest(m1.0)
bgtest(m1.1)
bgtest(m1.2)

acf(resid(m1.0),1000)
acf(resid(m1.1),1000)
acf(resid(m1.2),1000)

# AIC
AIC(m1.0,m1.1,m1.2)

# Correlation
cor(predict(m1.0,type="response"),y=Mon$MeanLepTotal,method='pearson')
# 0.4681963
cor(predict(m1.1,type="response"),y=Mon$MeanLepTotal,method='pearson')
# 0.6125846
cor(predict(m1.2,type="response"),y=Mon$MeanLepTotal,method='pearson')
# 0.510016

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results
summary(m1.0)

# Year and Time in the sea
newdata2 <- data.frame(
	Class=rep('S1',116),CageType=rep('Std',116),BayName=rep('Kilkieran',116),
	Year_f=factor(rep(1:29,each=4),levels=1:29,labels=levels(Mon$Year_f)),
	NewMonth=rep('Aug',116),
	TimeSea_f=factor(rep(1:4,29),levels = 1:4, labels =levels(Mon$TimeSea_f)))

newdata2 <- cbind(newdata2, predict(m1.0, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

newdata2$Year<-rep(1994:2022,each=4)

A<-ggplot(newdata2,aes(Year,Abun,colour=TimeSea_f)) +
	geom_point()+
	geom_line(,size = 1.1)+
	ggtitle('Pen level')+
	geom_ribbon(aes(ymin=LL,ymax=UL,fill=TimeSea_f), alpha = 0.3)+
	theme(text = element_text(size = 17))
A

# New Month, Time in sea
newdata2 <- data.frame(
	Class=rep('S1',116),CageType=rep('Std',116),BayName=rep('Kilkieran',116),
	Year_f=rep('2007',116),
	NewMonth=factor(rep(1:14,each=116),levels = 1:14,labels=levels(as.factor(Mon$NewMonth))),
	TimeSea_f=factor(rep(rep(1:4,29),14),levels=1:4,labels=levels(Mon$TimeSea_f)))

newdata2 <- cbind(newdata2, predict(m1.0, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2$NewMonth<-factor(newdata2$NewMonth,levels=c("Jan/Dec","Feb","Mar1",
	"Mar2","Apr1","Apr2","May1","May2","Jun","Jul","Aug","Sep","Oct","Nov"))

newdata2$Month<-rep(1:14,each=116)

B<-ggplot(newdata2,aes(Month,Abun,colour=TimeSea_f)) +
	geom_point()+
	geom_line(,size = 1.1)+
	ggtitle('Pen level')+
	geom_ribbon(aes(ymin=LL,ymax=UL,fill=TimeSea_f), alpha = 0.3)+
	theme(text = element_text(size = 17))
B

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################################################################################################################

############################# Farm level

# New dataset on the level of farm: mean per farm

Moi<-aggregate(Mon$MeanLepTotal,list(
	Mon$Year_n,Mon$Year_f,Mon$SiteName,Mon$NewMonth,Mon$Month,
	Mon$Class,Mon$TimeSea_f,Mon$SplDate,Mon$BayName),mean)
colnames(Moi)<-c('Year','Year_f','SiteName','NewMonth','Month','Class',
	'TimeSea','Date','BayName','Lep')
head(Moi)

# New variables
Moi$TimeSea<-as.factor(Moi$TimeSea)

# Exploration of autocorrelation for data: ACF
Moi$Date<-paste0(Moi$Year_f,'/',Moi$Month,'/01')
Moi$Date<-as.Date(Moi$Date)
Mmm<-xts(Moi$Lep,Moi$Date) 

acf(coredata(Mmm),1000)

# Full Model
m2.0<-glmmTMB(Lep~Class+Year_f+NewMonth+TimeSea+BayName+(1|SiteName),
		ziformula=~Class+Year_f+NewMonth+TimeSea+BayName,
		family=ziGamma(link="log"),data=Moi)
summary(m2.0)

# Checking autocorrelation and comparison of models
durbinWatsonTest(resid(m2.0))
dwtest(m2.0)
bgtest(m2.0)

acf(resid(m2.0),1000)

# Correlation
cor(predict(m2.0,type="response"),y=Moi$Lep,method='pearson')
# 0.466228

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results
summary(m2.0)

# Year
newdata2 <- data.frame(
	Class=rep('S1',116),BayName=rep('Kilkieran',116),
	Year_f=factor(rep(1:29,each=4),levels=1:29,labels=levels(Moi$Year_f)),
	NewMonth=rep('Aug',116),
	TimeSea=factor(rep(1:4,29),levels = 1:4, labels =levels(Moi$TimeSea)))

newdata2 <- cbind(newdata2, predict(m2.0, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2$Year<-rep(1994:2022,each=4)

C<-ggplot(newdata2,aes(Year,Abun,colour=TimeSea)) +
	geom_point()+
	geom_line(,size = 1.1)+
	ggtitle('Farm level')+
	geom_ribbon(aes(ymin=LL,ymax=UL,fill=TimeSea), alpha = 0.3)+
	theme(text = element_text(size = 17))
C

# New Month
newdata2 <- data.frame(
	Class=rep('S1',116),CageType=rep('Std',116),BayName=rep('Kilkieran',116),
	Year_f=rep('2007',116),
	NewMonth=factor(rep(1:14,each=116),levels = 1:14,labels=levels(as.factor(Moi$NewMonth))),
	TimeSea=factor(rep(rep(1:4,29),14),levels=1:4,labels=levels(Moi$TimeSea)))

newdata2 <- cbind(newdata2, predict(m2.0, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2$NewMonth<-factor(newdata2$NewMonth,levels=c("Jan/Dec","Feb","Mar1",
	"Mar2","Apr1","Apr2","May1","May2","Jun","Jul","Aug","Sep","Oct","Nov"))

newdata2$Month<-rep(1:14,each=116)

D<-ggplot(newdata2,aes(Month,Abun,colour=TimeSea)) +
	geom_point()+
	geom_line(,size = 1.1)+
	ggtitle('Farm level')+
	geom_ribbon(aes(ymin=LL,ymax=UL,fill=TimeSea), alpha = 0.3)+
	theme(text = element_text(size = 17))
D

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################################################################################################

# Comparison of variation between datasets

# Pen level
sd(Mon$MeanLepTotal)
mean(Mon$MeanLepTotal)
median(Mon$MeanLepTotal)
summary(Mon$MeanLepTotal)

# Farm level
sd(Moi$Lep)
mean(Moi$Lep)
median(Moi$Lep)
summary(Moi$Lep)

# Number of observation per dataset
# Pen level
nrow(Mon)
# Farm level
nrow(Moi)

# Comparison of models

# Pen level
summary(m1.0)

# Farm level
summary(m2.0)

# Comparison of results
grid.arrange(A,C,ncol = 2)
grid.arrange(B,D,ncol = 2)


