
# Task 1
# Question 2: The effect of the location of salmon farms on sea lice abundance

library(ggplot2)
library(glmmTMB)
library(xts)
library(car)
library(lmtest)
library(gridExtra)

# Sea lice database
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
Moneo<-Moneo[Moneo$Zdroj %in% 'Individ',]

# Adjustment of variables
Moneo$Year_f<-as.factor(Moneo$Year)
Moneo$TimeSea_f<-as.factor(Moneo$TimeSea)

Moneo[!(Moneo$New_Comp %in% c("Mowi","Bifand/Mowi","Bradan Beo Teo.")),]$New_Comp<-'Others'
Moneo[Moneo$New_Comp %in% c("Bifand/Mowi"),]$New_Comp<-"Mowi"

# Adjustment of Pen as a variable with more levels
AA<-aggregate(Moneo$S.R,list(Moneo$SiteName,Moneo$NewMonth,Moneo$Year,Moneo$S.R,Moneo$SampleID),length)

AA$Group.6<-paste(AA$Group.5,AA$Group.4)
AA$Group.7<-paste(AA$Group.1,AA$Group.2,AA$Group.3)
AA$Group.8<-NA
tapply(AA$Group.6,paste(AA$Group.1,AA$Group.2,AA$Group.3),unique)

for(i in 1:length(unique(AA$Group.7))){
	A<-AA[AA$Group.7 %in% unique(AA$Group.7)[i],]
	B<-A[A$Group.4 %in% c('Random'),]
	C<-A[A$Group.4 %in% c('Standard'),]
	if(nrow(B)>0){for (e in 1:nrow(B)){
		Id<-B[e,]$Group.5
		AA[AA$Group.5 %in% Id,]$Group.8<-paste0(B[e,]$Group.4,e)
	}}
	if(nrow(C)>0){for (j in 1:nrow(C)){
		Id<-C[j,]$Group.5
		AA[AA$Group.5 %in% Id,]$Group.8<-paste0(C[j,]$Group.4,j)
	}}
	}

# Combining datasets
colnames(AA)<-c('SiteName','NewMonth','Year','S.R','SampleID','Count','ID1','ID2','New_Pen1')
AA_F<-AA[,c(1,2,3,4,5,9)]
Monee<-merge(Moneo,AA_F,by=c('SiteName','NewMonth','Year','S.R','SampleID'))

# Discarding outliers
Mon_1<-Monee[Monee$LepeophtheirusSalmonisTotal<161,]

Mon_1[c(13527), ]
Mon_1[c(71131), ]
Mon_1.1<- Mon_1[-c(13527,71131), ]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################################################################################################################

############################# Fish level

# New variables
Mon_1.1$Year_n<-as.numeric(Mon_1.1$Year_f)/10
Mon_1.1$Year_w<-numFactor(Mon_1.1$Year)
Mon_1.1$Month_w<-numFactor(Mon_1.1$Month)
Mon_1.1$New_Pen1<-as.numeric(as.factor(Mon_1.1$New_Pen1))

# Exploration of autocorrelation for data
Mon_1.1$Date<-as.Date(Mon_1.1$Date)   
Mmm<-xts(Mon_1.1$LepeophtheirusSalmonisTotal,Mon_1.1$Date) 

acf(coredata(Mmm),1000)

# Without covariate structure for temporal autocorrelation
m1<-glmmTMB(LepeophtheirusSalmonisTotal~scale(Depth_FINAL)+
	scale(Distance_Mouth)+scale(Length)+scale(Orientation)+
	scale(Width_Mouth)+Year_n+NewMonth+TimeSea_f+
	Class+New_Comp+scale(Fallow)+(1|New_Pen1),
	family=nbinom2,data=Mon_1.1)
summary(m1)

# Covariate structure: Year_w
m1.1<-glmmTMB(LepeophtheirusSalmonisTotal~scale(Depth_FINAL)+
	scale(Distance_Mouth)+scale(Length)+scale(Orientation)+
	scale(Width_Mouth)+Year_n+NewMonth+TimeSea_f+
	Class+New_Comp+scale(Fallow)+(1|New_Pen1)+ou(Year_w+0|New_Pen1),
	family=nbinom2,data=Mon_1.1)
summary(m1.1)

# Covariate structure: Month_w
m1.2<-glmmTMB(LepeophtheirusSalmonisTotal~scale(Depth_FINAL)+
	scale(Distance_Mouth)+scale(Length)+scale(Orientation)+
	scale(Width_Mouth)+Year_n+NewMonth+TimeSea_f+
	Class+New_Comp+scale(Fallow)+(1|New_Pen1)+ou(Month_w+0|New_Pen1),
	family=nbinom2,data=Mon_1.1)
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

# Correlation
cor(predict(m1,type="response"),y=Mon_1.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3893499
cor(predict(m1.1,type="response"),y=Mon_1.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.426656
cor(predict(m1.2,type="response"),y=Mon_1.1$LepeophtheirusSalmonisTotal,method='pearson')
# 0.3961354

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results
# only variables of interest
summary(m1)

# Depth and Time in the sea
newdata2<- data.frame(
	Depth_FINAL= rep(seq(from = min(Mon_1.1$Depth_FINAL), to = max(Mon_1.1$Depth_FINAL),length.out = 100),4),
	Distance_Mouth= mean(scale(Mon_1.1$Distance_Mouth)),
	Length= mean(scale(Mon_1.1$Length)),Orientation= mean(scale(Mon_1.1$Orientation)),
	Width_Mouth= mean(scale(Mon_1.1$Width_Mouth)),Year_n= mean(as.vector(Mon_1.1$Year_n)),
	NewMonth= rep('Nov',400),	
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_1.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),Fallow= mean(scale(Mon_1.1$Fallow)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

A<-ggplot(newdata2, aes(Depth_FINAL,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Depth_FINAL", y = "Abundance")+ 
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Distance and Time in the sea
newdata2<- data.frame(
	Depth_FINAL= mean(scale(Mon_1.1$Depth_FINAL)),
	Distance_Mouth= rep(seq(from = min(Mon_1.1$Distance_Mouth), to = max(Mon_1.1$Distance_Mouth),length.out = 100),4),
	Length= mean(scale(Mon_1.1$Length)),Orientation= mean(scale(Mon_1.1$Orientation)),
	Width_Mouth= mean(scale(Mon_1.1$Width_Mouth)),Year_n= mean(as.vector(Mon_1.1$Year_n)),
	NewMonth= rep('Nov',400),	
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_1.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),Fallow= mean(scale(Mon_1.1$Fallow)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

B<-ggplot(newdata2, aes(Distance_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Distance_Mouth", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Length and Time in the sea
newdata2<- data.frame(
	Depth_FINAL= mean(scale(Mon_1.1$Depth_FINAL)),
	Distance_Mouth= mean(scale(Mon_1.1$Distance_Mouth)),
	Length= rep(seq(from = min(Mon_1.1$Length), to = max(Mon_1.1$Length),length.out = 100),4),
	Orientation= mean(scale(Mon_1.1$Orientation)),
	Width_Mouth= mean(scale(Mon_1.1$Width_Mouth)),Year_n= mean(as.vector(Mon_1.1$Year_n)),
	NewMonth= rep('Nov',400),	
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_1.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),Fallow= mean(scale(Mon_1.1$Fallow)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

C<-ggplot(newdata2, aes(Length,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Length", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Orientation and Time in the sea
newdata2<- data.frame(
	Depth_FINAL= mean(scale(Mon_1.1$Depth_FINAL)),
	Distance_Mouth= mean(scale(Mon_1.1$Distance_Mouth)),
	Length= mean(scale(Mon_1.1$Length)),
	Orientation= rep(seq(from = min(Mon_1.1$Orientation), to = max(Mon_1.1$Orientation),length.out = 100),4),
	Width_Mouth= mean(scale(Mon_1.1$Width_Mouth)),Year_n= mean(as.vector(Mon_1.1$Year_n)),
	NewMonth= rep('Nov',400),	
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_1.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),Fallow= mean(scale(Mon_1.1$Fallow)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

D<-ggplot(newdata2, aes(Orientation,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Orientation", y = "Abundance")+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Width_Mouth and Time in the sea
newdata2<- data.frame(
	Depth_FINAL= mean(scale(Mon_1.1$Depth_FINAL)),
	Distance_Mouth= mean(scale(Mon_1.1$Distance_Mouth)),
	Length= mean(scale(Mon_1.1$Length)),
	Orientation= mean(scale(Mon_1.1$Orientation)),
	Width_Mouth= rep(seq(from = min(Mon_1.1$Width_Mouth), to = max(Mon_1.1$Width_Mouth),length.out = 100),4),
	Year_n= mean(as.vector(Mon_1.1$Year_n)),
	NewMonth= rep('Nov',400),	
	TimeSea_f = factor(rep(1:4, each = 100), levels = 1:4, labels =levels(Mon_1.1$TimeSea_f)),
	Class = rep('S1',400),New_Comp = rep('Mowi',400),Fallow= mean(scale(Mon_1.1$Fallow)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

E<-ggplot(newdata2, aes(Width_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea_f), alpha = .25) +
  geom_line(aes(colour = TimeSea_f), linewidth = 2) +
  labs(x = "Width_Mouth", y = "Abundance")+ 
	ylab('')+
  ggtitle('Fish level')+
	theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(A,B,C,D,E,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

############################################################################################################################################

############################# Pen level

# New dataset on the level of pen: mean per pen

Mon<-aggregate(Mon_1.1$LepeophtheirusSalmonisTotal,list(
	Mon_1.1$Year_n,Mon_1.1$Year,Mon_1.1$SiteName,Mon_1.1$Month,
	Mon_1.1$NewMonth,Mon_1.1$New_Pen1,Mon_1.1$Class,Mon_1.1$TimeSea,
	Mon_1.1$New_Comp,Mon_1.1$SampleID,Mon_1.1$Depth_FINAL,Mon_1.1$Distance_Mouth,
	Mon_1.1$Length,Mon_1.1$Orientation,Mon_1.1$Width_Mouth,Mon_1.1$Fallow),mean)
colnames(Mon)<-c('Year','Year_f','SiteName','Month','NewMonth','New_Pen','Class','TimeSea',
	'Company','ID','Depth','Distance_Mouth','Length','Orientation','Width_Mouth','Fallow','Lep')
head(Mon)

# New variables
Mon$TimeSea<-as.factor(Mon$TimeSea)

# Exploration of autocorrelation for data
Mon$Date<-paste0(Mon$Year_f,'/',Mon$Month,'/01')
Mon$Date<-as.Date(Mon$Date)
Mmm<-xts(Mon$Lep,Mon$Date) 

acf(coredata(Mmm),1000)

# Full Model
m2<-glmmTMB(Lep~scale(Depth)+scale(Distance_Mouth)+scale(Length)+
	scale(Orientation)+scale(Width_Mouth)+Year+NewMonth+TimeSea+
	Class+Company+scale(Fallow)+(1|New_Pen),
	ziformula=~scale(Depth)+scale(Distance_Mouth)+scale(Length)+
	scale(Orientation)+scale(Width_Mouth)+Year+NewMonth+TimeSea+
	Class+Company+scale(Fallow),
	family=ziGamma(link="log"),data=Mon)
summary(m2)

# Checking autocorrelation and comparison of models
durbinWatsonTest(resid(m2))
dwtest(m2)
bgtest(m2)

acf(resid(m2),1000)

# Correlation
cor(predict(m2,type="response"),y=Mon$Lep,method='pearson')
# 0.4452394

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results
# only looking at the variables of interest
summary(m2)

# Depth and Time in the sea
newdata2<- data.frame(
	Depth= rep(seq(from = min(Mon$Depth),to=max(Mon$Depth),length.out=100),4),
	Distance_Mouth= mean(scale(Mon$Distance_Mouth)),
	Length= mean(scale(Mon$Length)),Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

F<-ggplot(newdata2, aes(Depth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Depth", y = "Abundance")+ 
  ggtitle('Pen level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Distance and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=rep(seq(from=min(Mon$Distance_Mouth),to=max(Mon$Distance_Mouth),length.out=100),4),
	Length= mean(scale(Mon$Length)),Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

G<-ggplot(newdata2, aes(Distance_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Distance_Mouth", y = "Abundance")+ 
  ggtitle('Pen level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

# Length and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length= rep(seq(from=min(Mon$Length),to=max(Mon$Length),length.out=100),4),
	Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

H<-ggplot(newdata2, aes(Length,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Length", y = "Abundance")+ 
  ggtitle('Pen level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

# Orientation and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length=mean(scale(Mon$Length)),
	Orientation= rep(seq(from=min(Mon$Orientation),to=max(Mon$Orientation),length.out=100),4),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

I<-ggplot(newdata2, aes(Orientation,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Orientation", y = "Abundance")+
  ggtitle('Pen level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Width_Mouth and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length=mean(scale(Mon$Length)),
	Orientation=mean(scale(Mon$Orientation)),
	Width_Mouth=rep(seq(from=min(Mon$Width_Mouth),to=max(Mon$Width_Mouth),length.out=100),4),
	Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

J<-ggplot(newdata2, aes(Width_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Width_Mouth", y = "Abundance")+ 
  ggtitle('Pen level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(F,G,H,I,J,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

############################################################################################################################################

############################# Farm level

# New dataset on the level of farm: mean per farm

Moi<-aggregate(Mon_1.1$LepeophtheirusSalmonisTotal,list(
	Mon_1.1$Year_n,Mon_1.1$Year,Mon_1.1$SiteName,Mon_1.1$Month,
	Mon_1.1$NewMonth,Mon_1.1$Class,Mon_1.1$TimeSea,Mon_1.1$Date,
	Mon_1.1$New_Comp,Mon_1.1$Depth_FINAL,Mon_1.1$Distance_Mouth,
	Mon_1.1$Length,Mon_1.1$Orientation,Mon_1.1$Width_Mouth,Mon_1.1$Fallow),mean)
colnames(Moi)<-c('Year','Year_f','SiteName','Month','NewMonth','Class','TimeSea',
	'Date','Company','Depth','Distance_Mouth','Length','Orientation',
	'Width_Mouth','Fallow','Lep')
head(Moi)

# New variables
Moi$TimeSea<-as.factor(Moi$TimeSea)

# Month weird format
# Exploration of autocorrelation for data
Moi$Date<-paste0(Moi$Year_f,'/',Moi$Month,'/01')
Moi$Date<-as.Date(Moi$Date)
Mmm<-xts(Moi$Lep,Moi$Date) 

acf(coredata(Mmm),1000)

# Full Model
m3<-glmmTMB(Lep~scale(Depth)+scale(Distance_Mouth)+scale(Length)+
	scale(Orientation)+scale(Width_Mouth)+Year+NewMonth+TimeSea+
	Class+Company+scale(Fallow),
	ziformula=~scale(Depth)+scale(Distance_Mouth)+scale(Length)+
	scale(Orientation)+scale(Width_Mouth)+Year+NewMonth+TimeSea+
	Class+Company+scale(Fallow),
	family=ziGamma(link="log"),data=Moi)
summary(m3)

# Checking autocorrelation and comparison of models
durbinWatsonTest(resid(m3))
dwtest(m3)
bgtest(m3)

acf(resid(m3),10000)

# Correlation
cor(predict(m3,type="response"),y=Moi$Lep,method='pearson')
# 0.4743219

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Exploration of results
summary(m3)

newdata1.1 <- data.frame(
	Depth= mean(scale(Mon$Depth)),Distance_Mouth= mean(scale(Mon$Distance_Mouth)),
	Length= mean(scale(Mon$Length)),Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= factor(rep(1:14,each=24), levels = 1:14, labels =levels(as.factor(Mon$NewMonth))),	
	TimeSea = factor(rep(rep(rep(1:4,times=2),each=3),14), levels = 1:4, labels =levels(Mon$TimeSea)),
	Class = factor(rep(rep(1:2, each = 12),14), levels = 1:2, labels =levels(as.factor(Mon$Class))),
	Company = factor(rep(rep(1:3,8),14), levels = 1:3, labels =levels(as.factor(Mon$Company))),
	Fallow= mean(scale(Mon$Fallow)))

newdata1.1$phat <- predict(m3,newdata1.1,type = "response",re.form=NA)
newdata1.1

# Depth and Time in the sea
newdata2<- data.frame(
	Depth= rep(seq(from = min(Mon$Depth),to=max(Mon$Depth),length.out=100),4),
	Distance_Mouth= mean(scale(Mon$Distance_Mouth)),
	Length= mean(scale(Mon$Length)),Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

K<-ggplot(newdata2, aes(Depth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Depth", y = "Abundance")+ 
  ggtitle('Farm level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Distance and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=rep(seq(from=min(Mon$Distance_Mouth),to=max(Mon$Distance_Mouth),length.out=100),4),
	Length= mean(scale(Mon$Length)),Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

L<-ggplot(newdata2, aes(Distance_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Distance_Mouth", y = "Abundance")+ 
  ggtitle('Farm level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

# Length and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length= rep(seq(from=min(Mon$Length),to=max(Mon$Length),length.out=100),4),
	Orientation= mean(scale(Mon$Orientation)),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

M<-ggplot(newdata2, aes(Length,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Length", y = "Abundance")+ 
  ggtitle('Farm level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

# Orientation and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length=mean(scale(Mon$Length)),
	Orientation= rep(seq(from=min(Mon$Orientation),to=max(Mon$Orientation),length.out=100),4),
	Width_Mouth= mean(scale(Mon$Width_Mouth)),Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

N<-ggplot(newdata2, aes(Orientation,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Orientation", y = "Abundance")+
  ggtitle('Farm level')+
	theme(legend.position = "none",text = element_text(size = 17))

# Width_Mouth and Time in the sea
newdata2<- data.frame(
	Depth=mean(scale(Mon$Depth)),
	Distance_Mouth=mean(scale(Mon$Distance_Mouth)),
	Length=mean(scale(Mon$Length)),
	Orientation=mean(scale(Mon$Orientation)),
	Width_Mouth=rep(seq(from=min(Mon$Width_Mouth),to=max(Mon$Width_Mouth),length.out=100),4),
	Year= mean(as.vector(Mon$Year)),
	NewMonth= rep('Nov',400),	
	TimeSea = factor(rep(1:4,each=100),levels=1:4,labels=levels(Mon$TimeSea)),
	Class = rep('S1',400),Company = rep('Mowi',400),Fallow= mean(scale(Mon$Fallow)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE,re.form=NA))
newdata2 <- within(newdata2, {
  Abun<- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

O<-ggplot(newdata2, aes(Width_Mouth,Abun)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = TimeSea), alpha = .25) +
  geom_line(aes(colour = TimeSea), linewidth = 2) +
  labs(x = "Width_Mouth", y = "Abundance")+ 
  ggtitle('Farm level')+
	ylab('')+
	theme(legend.position = "none",text = element_text(size = 17))

grid.arrange(K,L,M,N,O,ncol = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################################################################################################

# Comparison of variation between datasets

# Fish level
sd(Mon_1.1$LepeophtheirusSalmonisTotal)
mean(Mon_1.1$LepeophtheirusSalmonisTotal)
median(Mon_1.1$LepeophtheirusSalmonisTotal)
summary(Mon_1.1$LepeophtheirusSalmonisTotal)

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
nrow(Mon_1.1)
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

# Each variable separately
grid.arrange(A,F,K,ncol = 3)
grid.arrange(B,G,L,ncol = 3)
grid.arrange(C,H,M,ncol = 3)
grid.arrange(D,I,N,ncol = 3)
grid.arrange(E,J,O,ncol = 3)




