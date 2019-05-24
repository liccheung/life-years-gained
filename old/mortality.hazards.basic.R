rm(list=ls())
library(survival)
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis_imputed1.RData")

nhis <- subset(nhis,nhis$age>=40)
nhis$deathage <- nhis$deathyears+nhis$age
morat <- coxph(Surv(age, deathage, died) ~ female, data = nhis, model = TRUE, y = TRUE) 
basehaz_morat <- basehaz(morat)

lcd <- coxph(Surv(age, deathage, lung.cancer.death) ~ female, data = nhis, model = TRUE, y = TRUE) 
basehaz_lcd <- basehaz(lcd)


plot(basehaz_morat$time,basehaz_morat$hazard,type="l",main="Mortality hazards in NHIS 97-01",xlab="age",ylab="Cumulative Hazard")             
lines(basehaz_morat$time,basehaz_morat$hazard*exp(morat$coefficients),col="green")
lines(basehaz_lcd$time,basehaz_lcd$hazard,col="red")
lines(basehaz_lcd$time,basehaz_lcd$hazard*exp(lcd$coefficients),col="blue")
legend(40,2.5,c("Male overall deaths","Female overall deaths","Male lung cancer deaths","Female lung cancer deaths"),
      cex=0.8,col=c("black","green","red","blue"),lty=c(1,1,1,1))