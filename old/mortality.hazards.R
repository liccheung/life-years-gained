rm(list=ls())
library(survival)
library(survey)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")

nhis1 <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis1, nest=TRUE)
design <- subset(master, age>=40)
morat <- svycoxph(Surv(age, deathage, died) ~ female+current, design=design, data = nhis1, model = TRUE, y = TRUE) 
basehaz_morat <- basehaz(morat)
basehaz_morat <- basehaz_morat[(basehaz_morat$time-floor(basehaz_morat$time))==.125,]
basehaz_morat$diff <- c(0,diff(basehaz_morat$hazard))
lm(formula=diff~time,basehaz_morat)  #0.005617


nhis2 <- subset(nhis,deathage2>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis2, nest=TRUE)
design <- subset(master, age>=40)
lcd <- svycoxph(Surv(age, deathage2, lung.cancer.death) ~ female+current, design=design, data = nhis2, model = TRUE, y = TRUE) 
basehaz_lcd <- basehaz(lcd)
basehaz_lcd <- basehaz_lcd[(basehaz_lcd$time-floor(basehaz_lcd$time))==.125,]
basehaz_lcd$diff <- c(0,diff(basehaz_lcd$hazard))
lm(formula=diff~time,basehaz_lcd)  #0.0001227

plot(basehaz_morat$time,basehaz_morat$diff,type="l",main="Mortality hazards in NHIS 97-01",xlab="age",ylab="Hazard")             
lines(basehaz_morat$time,basehaz_morat$diff*exp(morat$coefficients[1]),col="green")
lines(basehaz_morat$time,basehaz_morat$diff*exp(morat$coefficients[2]),col="red")
lines(basehaz_morat$time,basehaz_morat$diff*exp(morat$coefficients[1]+morat$coefficients[2]),col="blue")
lines(basehaz_lcd$time,basehaz_lcd$diff,type="l",lty=2,col="black")
lines(basehaz_lcd$time,basehaz_lcd$diff*exp(morat$coefficients[1]),type="l",lty=2,col="green")
lines(basehaz_lcd$time,basehaz_lcd$diff*exp(morat$coefficients[2]),type="l",lty=2,col="red")
lines(basehaz_lcd$time,basehaz_lcd$diff*exp(morat$coefficients[1]+morat$coefficients[2]),type="l",lty=2,col="blue")
legend(40,3.5,c("Male overall deaths","Female overall deaths","Male lung cancer deaths","Female lung cancer deaths"),
      cex=0.8,col=c("black","green","red","blue"),lty=c(1,1,1,1))