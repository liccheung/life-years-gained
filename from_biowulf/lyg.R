library(survival)
library(survey)
load(file="/home/cheunglc/lyg/mortality.model.v4.RData")
source("/home/cheunglc/lyg/lcrisk.R")

lyg <- function(x){
  print(x)
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  ely <- sum(diff(c(z$age[x],c$time))*c$surv)  
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x]+5,stop.time=100)
  lyg <- (1-z$lcd1[x])**.8-(1-z$lcd1[x]) + 
    (1-z$lcd2[x])**.8-(1-z$lcd2[x]) +
    (1-z$lcd3[x])**.8-(1-z$lcd3[x]) +
    (1-z$lcd4[x])**.8-(1-z$lcd4[x]) +
    ((1-z$lcd5[x])**.8-(1-z$lcd5[x]))*(1+sum(diff(c(z$age[x]+5,c$time))*c$surv))
  res <- c(ely,ely+lyg,lyg)
  return(res)
}

persons <- data.frame(age=y$age,
                      female=y$female,
                      smkyears=y$smkyears,
                      qtyears=y$qtyears,
                      cpd=y$cpd,
                      race=as.numeric(as.character(y$race)),
                      emp=y$emp,
                      fam.lung.trend=y$fam.lung.trend,
                      bmi=y$bmi,
                      edu6=y$edu6)

missval <- ifelse(!is.na(persons$age) & !is.na(persons$female) & !is.na(persons$smkyears) & !is.na(persons$qtyears) &
       	       	  !is.na(persons$cpd) & !is.na(persons$race) & !is.na(persons$emp) & 
                  !is.na(persons$fam.lung.trend) & !is.na(persons$bmi) & !is.na(persons$edu6),0,1)

y$lcd1[missval==0] <- lcrisk(persons[missval==0,],1)[,3]/1000
y$lcd2[missval==0] <- lcrisk(persons[missval==0,],2)[,3]/1000   
y$lcd3[missval==0] <- lcrisk(persons[missval==0,],3)[,3]/1000   
y$lcd4[missval==0] <- lcrisk(persons[missval==0,],4)[,3]/1000   
y$lcd5[missval==0] <- lcrisk(persons[missval==0,],5)[,3]/1000   

y$birthyear <- y$year-y$age

z <- subset(y,missval==0)
lyg_res <- sapply(1:nrow(z),lyg)
lyg_res <- t(lyg_res)
y$ely[missval==0] <- lyg_res[,1]
y$ely_ct[missval==0] <- lyg_res[,2]
y$lyg[missval==0] <- lyg_res[,3]
