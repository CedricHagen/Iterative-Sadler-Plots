#This script generates the plots discussed in
#Hagen (in review, The Sedimentary Record). There
#is also a function so that these small analyses 
#can be repeated for other age-depth models.

#Last Modified on July 3, 2024

#Author: Dr. Cedric Hagen, University of Colorado, Boulder  

# Synthetic 1 Setup ------------------------------------------------------

age1 <- 530
age2 <- 525
age3 <- 519
age4 <- 518

height1 <- 0
height2 <- 33
height3 <- 67
height4 <- 100

test_height <- 1:99

test_age <- numeric(0)
age_range <- age1-age2
height_range <- height2-height1
heights <- 1:33
ages <- age1 -(heights/height_range)*age_range
test_age <- c(test_age,ages)

age_range <- age2-age3
height_range <- height3-height2
heights <- 1:33
ages <- age2 - (heights/height_range)*age_range
test_age <- c(test_age,ages)

age_range <- age3-age4
height_range <- height4-height3
heights <- 1:33
ages <- age3 - (heights/height_range)*age_range
test_age <- c(test_age,ages)

#plot(test_height,test_age,pch=20,ylim=rev(range(test_age)),ylab='Age (Ma)',xlab='Height (m)')
#lines(test_height,test_age, lwd=2)
#abline(h = age1, col="red", lty=2)
#abline(h = age2, col="red", lty=2)
#abline(h = age3, col="red", lty=2)
#abline(h = age4, col="red", lty=2)

# Synthetic 1: Plot 1 ------------------------------------------------------

durs <- numeric(0)
accums <- numeric(0)
iss <- numeric(0)
jss <- numeric(0)
oage <- numeric(0)
for (i in 1:length(test_age)) {
  for (j in 1:length(test_age)){
    if (j>1) {
      durs <- c(durs, test_age[i] - test_age[j])
      accums <- c(accums, (test_height[j] - test_height[i]) / (test_age[i] - test_age[j]))
      oage <- c(oage,test_age[j])
      iss <- c(iss, i)
      jss <- c(jss, j)
    }
  }
}

durs <- durs * 1000000  # yr
accums <- (accums * 1000) / 1000000  # mm/yr
inds <- which(durs>0)
durs <- durs[inds]
accums <- accums[inds]
oage <- oage[inds]

df = data.frame(durs,accums,oage)

library(ggplot2)
library(scales)

ggplot(df, aes(x=durs, y=accums, color=oage)) + 
  geom_point(size=2,pch=20) +
  scale_x_continuous(name = 'Duration (yr)', labels = label_log(digits = 2), limits = c(10^2, 10^8), trans='log10') +
  scale_y_continuous(name = 'Accumulation Rate (mm/yr)', labels = label_log(digits = 2), limits = c(10^-5, 10^2), trans='log10') +
  scale_color_viridis_c(option = "plasma", name = "Overlying\nAge (Ma)") +
  theme_classic()

# Synthetic 1: Plot 2 ------------------------------------------------------

library(drc)
library(nlme)
library(aomisc)
fit_result <- drm(accums ~ durs, fct = DRC.powerCurve())
summary(fit_result)

a <- fit_result$fit$par[1]
b <- fit_result$fit$par[2]

xval <- seq(min(durs),max(durs),10000)
yval <- a*(xval)^(b)

plot(fit_result, log="", main = "",xlab = "Duration (yr)", ylab = "Accumulation Rate (mm/yr)",pch=20)
lines(xval,yval,col="red",lwd=1.5,lty=1)
text(4e6, 0.04, paste("a = ",as.character(round(a,3)),", b = ",as.character(round(b,3))))

# Synthetic 1: Plot 3 ------------------------------------------------------

resolution <- seq(0.001, age1 - age4, 0.001)
complete <- (resolution / (age1 - age4)) ^ (b*-1)
ind <- which.min(abs(resolution - 1))
ind2 <- which.min(abs(complete - 0.50))
onemyrcomp <- complete[ind] * 100
fifty <- resolution[ind2]

plot(resolution, complete * 100, type = "l", col = "black", lwd = 2, xlab = "Completion interval (Myr)", ylab = "Stratigraphic completeness (%)",xlim = rev(range(resolution)))
abline(h = onemyrcomp, col = "red", lty = 2, lwd = 1.5)
abline(v = fifty, col = "red", lty = 2, lwd = 1.5)
axis(1, at = seq(0, age1 - age4, by = 10), labels = seq(0, age1 - age4, by = 10))
axis(2, at = seq(0, 100, by = 10))

# Synthetic 2 Setup ------------------------------------------------------

age1 <- 530
age2 <- 525
age3 <- 519
age4 <- 518

height1 <- 0
height2 <- 33
height3 <- 67
height4 <- 100

library(varbvs)
test_height <- 1
for (i in 1:98) {
  test_height <- c(test_height,test_height[end(test_height)[1]]+abs(randn(1,1)));
}

test_age <- numeric(0)
age_range <- age1-age2
height_range <- height2-height1
heights <- 1:33
ages <- age1 -(heights/height_range)*age_range
test_age <- c(test_age,ages)

age_range <- age2-age3
height_range <- height3-height2
heights <- 1:33
ages <- age2 - (heights/height_range)*age_range
test_age <- c(test_age,ages)

age_range <- age3-age4
height_range <- height4-height3
heights <- 1:33
ages <- age3 - (heights/height_range)*age_range
test_age <- c(test_age,ages)

plot(test_height,test_age,pch=20,ylim=rev(range(test_age)),ylab='Age (Ma)',xlab='Height (m)')
lines(test_height,test_age, lwd=2)
abline(h = age1, col="red", lty=2)
abline(h = age2, col="red", lty=2)
abline(h = age3, col="red", lty=2)
abline(h = age4, col="red", lty=2)

# Synthetic 2: Plot 1 ------------------------------------------------------

durs <- numeric(0)
accums <- numeric(0)
iss <- numeric(0)
jss <- numeric(0)
oage <- numeric(0)
for (i in 1:length(test_age)) {
  for (j in 1:length(test_age)){
    if (j>1) {
      durs <- c(durs, test_age[i] - test_age[j])
      accums <- c(accums, (test_height[j] - test_height[i]) / (test_age[i] - test_age[j]))
      oage <- c(oage,test_age[j])
      iss <- c(iss, i)
      jss <- c(jss, j)
    }
  }
}

durs <- durs * 1000000  # yr
accums <- (accums * 1000) / 1000000  # mm/yr
inds <- which(durs>0)
durs <- durs[inds]
accums <- accums[inds]
oage <- oage[inds]

df = data.frame(durs,accums,oage)

library(ggplot2)
library(scales)

ggplot(df, aes(x=durs, y=accums, color=oage)) + 
  geom_point(size=2,pch=20) +
  scale_x_continuous(name = 'Duration (yr)', labels = label_log(digits = 2), limits = c(10^2, 10^8), trans='log10') +
  scale_y_continuous(name = 'Accumulation Rate (mm/yr)', labels = label_log(digits = 2), limits = c(10^-5, 10^2), trans='log10') +
  scale_color_viridis_c(option = "plasma", name = "Overlying\nAge (Ma)") +
  theme_classic()

# Synthetic 2: Plot 2 ------------------------------------------------------

library(drc)
library(nlme)
library(aomisc)
fit_result <- drm(accums ~ durs, fct = DRC.powerCurve())
summary(fit_result)

a <- fit_result$fit$par[1]
b <- fit_result$fit$par[2]

xval <- seq(min(durs),max(durs),10000)
yval <- a*(xval)^(b)

plot(fit_result, log="", main = "",xlab = "Duration (yr)", ylab = "Accumulation Rate (mm/yr)",pch=20)
lines(xval,yval,col="red",lwd=1.5,lty=1)
text(4e6, 0.04, paste("a = ",as.character(round(a,3)),", b = ",as.character(round(b,3))))

# Synthetic 2: Plot 3 ------------------------------------------------------

resolution <- seq(0.001, age1 - age4, 0.001)
complete <- (resolution / (age1 - age4)) ^ (b*-1)
ind <- which.min(abs(resolution - 1))
ind2 <- which.min(abs(complete - 0.50))
onemyrcomp <- complete[ind] * 100
fifty <- resolution[ind2]

plot(resolution, complete * 100, type = "l", col = "black", lwd = 2, xlab = "Completion interval (Myr)", ylab = "Stratigraphic completeness (%)",xlim = rev(range(resolution)))
abline(h = onemyrcomp, col = "red", lty = 2, lwd = 1.5)
abline(v = fifty, col = "red", lty = 2, lwd = 1.5)
axis(1, at = seq(0, age1 - age4, by = 10), labels = seq(0, age1 - age4, by = 10))
axis(2, at = seq(0, 100, by = 10))


# Synthetic 3 Setup ------------------------------------------------------

test_height <- 1:99
test_age <- 530-2.5*log10(test_height)

plot(test_height,test_age,pch=20,ylim=rev(range(test_age)),ylab='Age (Ma)',xlab='Height (m)')
lines(test_height,test_age, lwd=2)
abline(h = 530, col="red", lty=2)

# Synthetic 3: Plot 1 ------------------------------------------------------

durs <- numeric(0)
accums <- numeric(0)
iss <- numeric(0)
jss <- numeric(0)
oage <- numeric(0)
for (i in 1:length(test_age)) {
  for (j in 1:length(test_age)){
    if (j>1) {
      durs <- c(durs, test_age[i] - test_age[j])
      accums <- c(accums, (test_height[j] - test_height[i]) / (test_age[i] - test_age[j]))
      oage <- c(oage,test_age[j])
      iss <- c(iss, i)
      jss <- c(jss, j)
    }
  }
}

durs <- durs * 1000000  # yr
accums <- (accums * 1000) / 1000000  # mm/yr
inds <- which(durs>0)
durs <- durs[inds]
accums <- accums[inds]
oage <- oage[inds]

df = data.frame(durs,accums,oage)

library(ggplot2)
library(scales)

ggplot(df, aes(x=durs, y=accums, color=oage)) + 
  geom_point(size=2,pch=20) +
  scale_x_continuous(name = 'Duration (yr)', labels = label_log(digits = 2), limits = c(10^2, 10^8), trans='log10') +
  scale_y_continuous(name = 'Accumulation Rate (mm/yr)', labels = label_log(digits = 2), limits = c(10^-5, 10^2), trans='log10') +
  scale_color_viridis_c(option = "plasma", name = "Overlying\nAge (Ma)") +
  theme_classic()

# Synthetic 3: Plot 2 ------------------------------------------------------

library(drc)
library(nlme)
library(aomisc)
fit_result <- drm(accums ~ durs, fct = DRC.powerCurve())
summary(fit_result)

a <- fit_result$fit$par[1]
b <- fit_result$fit$par[2]

xval <- seq(min(durs),max(durs),10000)
yval <- a*(xval)^(b)

plot(fit_result, log="", main = "",xlab = "Duration (yr)", ylab = "Accumulation Rate (mm/yr)",pch=20)
lines(xval,yval,col="red",lwd=1.5,lty=1)
text(4e6, 0.04, paste("a = ",as.character(round(a,3)),", b = ",as.character(round(b,3))))

# Synthetic 3: Plot 3 ------------------------------------------------------

resolution <- seq(0.001, age1 - age4, 0.001)
complete <- (resolution / (age1 - age4)) ^ (b*-1)
ind <- which.min(abs(resolution - 1))
ind2 <- which.min(abs(complete - 0.50))
onemyrcomp <- complete[ind] * 100
fifty <- resolution[ind2]

plot(resolution, complete * 100, type = "l", col = "black", lwd = 2, xlab = "Completion interval (Myr)", ylab = "Stratigraphic completeness (%)",xlim = rev(range(resolution)))
abline(h = onemyrcomp, col = "red", lty = 2, lwd = 1.5)
abline(v = fifty, col = "red", lty = 2, lwd = 1.5)
axis(1, at = seq(0, age1 - age4, by = 10), labels = seq(0, age1 - age4, by = 10))
axis(2, at = seq(0, 100, by = 10))


# Iterative Sadler Plot ------------------------------------------------------

durs <- numeric(0)
accums <- numeric(0)
iss <- numeric(0)
jss <- numeric(0)
oage <- numeric(0)
mean_sr<-numeric(0)
sd_sr<-numeric(0)
age_sr<-numeric(0)
mean_dur<-numeric(0)
sd_dur<-numeric(0)
for (i in 1:length(test_age)) {
  #print(i)
  sr<-numeric(0)
  dd<-numeric(0)
  for (j in i+1:length(test_age)){
    durs <- c(durs, test_age[i] - test_age[j])
    dd <- c(dd, test_age[i] - test_age[j])
    accums <- c(accums, (test_height[j] - test_height[i]) / (test_age[i] - test_age[j]))
    oage <- c(oage,test_age[j])
    iss <- c(iss, i)
    jss <- c(jss, j)
    sr <- c(sr,(test_height[j] - test_height[i]) / (test_age[i] - test_age[j]))
  }
  mean_sr<-c(mean_sr,mean(sr,na.rm=TRUE))
  sd_sr<-c(sd_sr,sd(sr,na.rm=TRUE))
  age_sr<-c(age_sr,test_age[i])
  mean_dur<-c(mean_dur,mean(dd,na.rm=TRUE))
  sd_dur<-c(sd_dur,sd(dd,na.rm=TRUE))
}


mean_sr <- (mean_sr * 1000) / 1000000  # mm/yr
sd_sr <- (sd_sr * 1000) / 1000000  # mm/yr
mean_dur <- mean_dur * 1000000  # yr
sd_dur <- sd_dur * 1000000  # yr

inds <- is.nan(mean_sr)
mean_sr[inds]<-NA
inds <- is.nan(sd_sr)
sd_sr[inds]<-NA

inds <- is.nan(mean_dur)
mean_dur[inds]<-NA
inds <- is.nan(sd_dur)
sd_dur[inds]<-NA

df = data.frame(mean_sr,sd_sr,mean_dur,sd_dur,age_sr)

library(ggplot2)
library(scales)

ggplot(df, aes(x=mean_dur, xmin=mean_dur-sd_dur, xmax=mean_dur+sd_dur, y=mean_sr, ymin=mean_sr-sd_sr, ymax=mean_sr+sd_sr, color=age_sr)) + 
  geom_point(size=2,pch=20) +
  geom_errorbar(aes(ymin=mean_sr-sd_sr, ymax=mean_sr+sd_sr),color='grey') +
  geom_errorbarh(aes(xmin=mean_dur-sd_dur, xmax=mean_dur+sd_dur),color='grey') +
  geom_point(size=2,pch=20) +
  scale_x_continuous(name = 'Duration (yr)', labels = label_log(digits = 2), limits = c(10^2, 10^8), trans='log10') +
  scale_y_continuous(name = 'Accumulation Rate (mm/yr)', labels = label_log(digits = 2), limits = c(10^-5, 10^2), trans='log10') +
  scale_color_viridis_c(option = "plasma", name = "Overlying\nAge (Ma)") +
  theme_classic()

# Long-term accumulation rate ------------------------------------------------------

df = data.frame(mean_sr,sd_sr,age_sr)

ggplot(df, aes(x=mean_sr, xmin=mean_sr-sd_sr, xmax=mean_sr+sd_sr, y=age_sr)) + 
  geom_point(size=2,pch=20,col='black') +
  geom_ribbon(alpha=0.5) +
  scale_x_continuous(name = 'Accumulation Rate (mm/yr)', limits = c(0, 1)) +
  scale_y_reverse(name = 'Age (Ma)', limits = c(max(test_age), min(test_age))) +
  theme_classic()

# Short-term accumulation rate ------------------------------------------------------


sr=numeric()
a=numeric()
for (i in 1:length(test_age)){
  if (i>1){
    sr <- c(sr, (test_height[i]-test_height[i-1])/(test_age[i-1] - test_age[i]))
    a <- c(a,test_age[i])
  }
}

sr <- (sr * 1000) / 1000000  # mm/yr

df = data.frame(sr,a,m = test_height[2:length(test_height)])

ggplot(df, aes(x=sr, y=a, color=m)) + 
  geom_point(size=2,pch=20) +
  scale_x_continuous(name = 'Accumulation Rate (mm/yr)', limits = c(0, 1)) +
  scale_y_reverse(name = 'Age (Ma)', limits = c(max(test_age), min(test_age))) +
  scale_color_viridis_c(option = "plasma", name = "Meterage") +
  theme_classic()



#This function will generate an 'iterative' Sadler plot for any age-depth model
#The function takes two variables: age, which is the age of each bed in Ma; and
#height, which is the stratigraphic height of each bed in meters. 
iterativeSadler <- function(age,height){
  
  durs <- numeric(0)
  accums <- numeric(0)
  iss <- numeric(0)
  jss <- numeric(0)
  oage <- numeric(0)
  mean_sr<-numeric(0)
  sd_sr<-numeric(0)
  age_sr<-numeric(0)
  mean_dur<-numeric(0)
  sd_dur<-numeric(0)
  for (i in 1:length(age)) {
    #print(i)
    sr<-numeric(0)
    dd<-numeric(0)
    for (j in i+1:length(age)){
      durs <- c(durs, age[i] - age[j])
      dd <- c(dd, age[i] - age[j])
      accums <- c(accums, (height[j] - height[i]) / (age[i] - age[j]))
      oage <- c(oage,age[j])
      iss <- c(iss, i)
      jss <- c(jss, j)
      sr <- c(sr,(height[j] - height[i]) / (age[i] - age[j]))
    }
    mean_sr<-c(mean_sr,mean(sr,na.rm=TRUE))
    sd_sr<-c(sd_sr,sd(sr,na.rm=TRUE))
    age_sr<-c(age_sr,age[i])
    mean_dur<-c(mean_dur,mean(dd,na.rm=TRUE))
    sd_dur<-c(sd_dur,sd(dd,na.rm=TRUE))
  }
  
  
  mean_sr <- (mean_sr * 1000) / 1000000  # mm/yr
  sd_sr <- (sd_sr * 1000) / 1000000  # mm/yr
  mean_dur <- mean_dur * 1000000  # yr
  sd_dur <- sd_dur * 1000000  # yr
  
  inds <- is.nan(mean_sr)
  mean_sr[inds]<-NA
  inds <- is.nan(sd_sr)
  sd_sr[inds]<-NA
  
  inds <- is.nan(mean_dur)
  mean_dur[inds]<-NA
  inds <- is.nan(sd_dur)
  sd_dur[inds]<-NA
  
  df = data.frame(mean_sr,sd_sr,mean_dur,sd_dur,age_sr)
  
  library(ggplot2)
  library(scales)
  
  ggplot(df, aes(x=mean_dur, xmin=mean_dur-sd_dur, xmax=mean_dur+sd_dur, y=mean_sr, ymin=mean_sr-sd_sr, ymax=mean_sr+sd_sr, color=age_sr)) + 
    geom_point(size=2,pch=20) +
    geom_errorbar(aes(ymin=mean_sr-sd_sr, ymax=mean_sr+sd_sr),color='grey') +
    geom_errorbarh(aes(xmin=mean_dur-sd_dur, xmax=mean_dur+sd_dur),color='grey') +
    geom_point(size=2,pch=20) +
    scale_x_continuous(name = 'Duration (yr)', labels = label_log(digits = 2), limits = c(10^2, 10^8), trans='log10') +
    scale_y_continuous(name = 'Accumulation Rate (mm/yr)', labels = label_log(digits = 2), limits = c(10^-5, 10^2), trans='log10') +
    scale_color_viridis_c(option = "plasma", name = "Overlying\nAge (Ma)") +
    theme_classic()
  
}

iterativeSadler(test_age,test_height)