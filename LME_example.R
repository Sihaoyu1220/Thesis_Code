# Load packages
library(nlme)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Load dataset
before_ART <- read.csv("before_ART.csv")
before_ART$PATIENT<-as.factor(before_ART$PATIENT)

# Keep first four observations for each patient
before_ART <- before_ART %>% 
  group_by(PATIENT) %>% 
  filter(row_number()<=4)
before_ART <- before_ART %>% group_by(PATIENT) %>% mutate(n=n())

# Randomly select four patients.
before_ART <- before_ART %>%
  filter(PATIENT==1|PATIENT==28|PATIENT==60|PATIENT==75)

# Fit LME model, with time as covariates. The LME has random intercept and slope.
long.dat <- groupedData ( log10 ~ time | PATIENT , data = before_ART )
lme <- lme(log10 ~ time, data = long.dat,random = ~ time|PATIENT,method="ML")
before_ART$predlm = predict(lme)

# Find fixed and random effects for each patient.
fixed <- summary(lme)$coefficients$fixed
ranef <- rownames_to_column(ranef(lme),"PATIENT")
subject1=before_ART[which(before_ART$PATIENT==2),]
subject2=before_ART[which(before_ART$PATIENT==46),]
subject3=before_ART[which(before_ART$PATIENT==60),]
subject4=before_ART[which(before_ART$PATIENT==75),]

# Make plot.
attach(mtcars)
par(mfrow=c(2,2))
attach(subject1)
sub1_value <- ranef[which(ranef$PATIENT==2),]
intercept <- as.numeric(sub1_value[2]+fixed[1])
slope <- as.numeric(sub1_value[3]+fixed[2])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylim=c(-0.5,7.5),ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(intercept+slope*x,from=0,to=28,lty=2,add=TRUE,lwd=2)
detach(subject1)
attach(subject2)
sub2_value <- ranef[which(ranef$PATIENT==46),]
intercept <- as.numeric(sub2_value[2]+fixed[1])
slope <- as.numeric(sub2_value[3]+fixed[2])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(intercept+slope*x,from=1,to=28,lty=2,add=TRUE,lwd=2)
detach(subject2)
attach(subject3)
sub3_value <- ranef[which(ranef$PATIENT==60),]
intercept <- as.numeric(sub3_value[2]+fixed[1])
slope <- as.numeric(sub3_value[3]+fixed[2])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylim=c(0,5),ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(intercept+slope*x,from=1,to=28,lty=2,add=TRUE,lwd=2)
detach(subject3)
attach(subject4)
sub4_value <- ranef[which(ranef$PATIENT==75),]
intercept <- as.numeric(sub4_value[2]+fixed[1])
slope <- as.numeric(sub4_value[3]+fixed[2])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylim=c(0,6),ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(intercept+slope*x,from=1,to=28,lty=2,add=TRUE,lwd=2)
detach(subject4)


