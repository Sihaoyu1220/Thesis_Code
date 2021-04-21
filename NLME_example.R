# Load package
library(nlme)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Load dataset.
before_ART <- read.csv("before_ART.csv")

# If the viral load values are censored, then add some noise on it.
data <- before_ART[! apply (is.na(before_ART ) ,1, any ) ,] #
n <- data %>% 
  filter(censor==1) %>% 
  nrow()
data <- data %>% 
  mutate(log10 =ifelse(censor==1, log10/2+rnorm(n,0,0.5), log10)) 
n2 <- data %>% 
  filter(log10<0) %>% 
  nrow()
while (n2>0){
data <- data %>% 
  mutate(log10 = ifelse(log10<0,0.5+rnorm(n2,0,0.25),log10))
n2 <- data %>% 
  filter(log10<0) %>% 
  nrow()
}
data <- data %>% 
  filter(PATIENT<50)
before_ART$PATIENT<-as.factor(before_ART$PATIENT)

# Fit NLME model.
long.dat <- groupedData ( log10 ~ time | PATIENT , data = before_ART )
logexp2 <- function (p1 ,b1 ,p2 ,b2 ,t) log10(exp(p1 -b1*t)+ exp (p2 -b2*t))
start <- c (p1=17 ,b1=4 ,p2=2.6 ,b2=0.05)
nlme.fit <- nlme(log10 ~ logexp2 (p1 ,b1 ,p2 ,b2 , time),
                     fixed = p1+b1+p2+b2 ~1,
                     random = p1+b1+p2+b2 ~1,
                     data = long.dat , start =c(start))
summary(nlme.fit)
before_ART$prednlm = predict(nlme.fit)
colnames(before_ART)[8]<-"predlm"

# Make plots.
ggplot(select, aes(x = time, y = log10, color = PATIENT) ) +
  geom_point() +
  geom_line(aes(y = predlm),size = 1)+
  theme_bw()+
  theme(text = element_text(size=14),legend.title=element_text(),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))

fixed <- summary(nlme.fit)$coefficients$fixed
ranef <- rownames_to_column(ranef(nlme.fit),"PATIENT")
subject1=before_ART[which(before_ART$PATIENT==5),]
subject2=before_ART[which(before_ART$PATIENT==9),]
subject3=before_ART[which(before_ART$PATIENT==60),]
subject4=before_ART[which(before_ART$PATIENT==75),]
attach(mtcars)
par(mfrow=c(2,2))
attach(subject1)
sub1_value <- ranef[which(ranef$PATIENT==5),]
p1 <- as.numeric(sub1_value[2]+fixed[1])
b1 <- as.numeric(sub1_value[3]+fixed[2])
p2 <- as.numeric(sub1_value[4]+fixed[3])
b2 <- as.numeric(sub1_value[5]+fixed[4])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(logexp2(p1,b1,p2,b2,x),from=1,to=26,lty=2,add=TRUE,lwd=2)
detach(subject1)
attach(subject2)
sub2_value <- ranef[which(ranef$PATIENT==9),]
p1 <- as.numeric(sub2_value[2]+fixed[1])
b1 <- as.numeric(sub2_value[3]+fixed[2])
p2 <- as.numeric(sub2_value[4]+fixed[3])
b2 <- as.numeric(sub2_value[5]+fixed[4])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(logexp2(p1,b1,p2,b2,x),from=1,to=21,lty=2,add=TRUE,lwd=2)
detach(subject2)
attach(subject3)
sub3_value <- ranef[which(ranef$PATIENT==60),]
p1 <- as.numeric(sub3_value[2]+fixed[1])
b1 <- as.numeric(sub3_value[3]+fixed[2])
p2 <- as.numeric(sub3_value[4]+fixed[3])
b2 <- as.numeric(sub3_value[5]+fixed[4])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(logexp2(p1,b1,p2,b2,x),from=1,to=21,lty=2,add=TRUE,lwd=2)
detach(subject3)
attach(subject4)
sub4_value <- ranef[which(ranef$PATIENT==75),]
p1 <- as.numeric(sub4_value[2]+fixed[1])
b1 <- as.numeric(sub4_value[3]+fixed[2])
p2 <- as.numeric(sub4_value[4]+fixed[3])
b2 <- as.numeric(sub4_value[5]+fixed[4])
plot(time,log10,pch=20,main="fitted vs. observed",xlab="time (in months)",ylab=bquote("Viral load (in" ~ log[10]~"-scale)"),mgp=c(2.5,1,0),cex=1.5,cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
curve(logexp2(p1,b1,p2,b2,x),from=1,to=25,lty=2,add=TRUE,lwd=2)
detach(subject4)
