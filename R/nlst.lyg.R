rm(list=ls(all=TRUE))
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/nlst/anlst_010918.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v6.RData")

lyg <- function(x){
  print(x)
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  ely <- sum(diff(c(z$age[x],c$time))*c$surv)
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x]+5,stop.time=100)
  lyg <- z$lcds1[x] + z$lcds2[x] + z$lcds3[x] + z$lcds4[x] +
    z$lcds5[x]*(1+sum(diff(c(z$age[x]+5,c$time))*c$surv))
  res <- c(ely,ely+lyg,lyg)
  return(res)
}

anlst$edu6 <- anlst$ed
anlst$year <- 2004
anlst$chd <- 0
anlst$angina <- 0
anlst$heartattack <- anlst$heart
anlst$heartdisease <- 0
anlst$kidney <- 0
anlst$liver <- 0
anlst$speceq <- 0
persons <- data.frame(age=anlst$age,
                      female=anlst$female,
                      smkyears=anlst$smkyears,
                      qtyears=anlst$qtyears,
                      cpd=anlst$cpd,
                      race=as.numeric(as.character(anlst$race)),
                      emp=anlst$emp,
                      fam.lung.trend=anlst$fam.lung.trend,
                      bmi=anlst$bmi,
                      edu6=anlst$edu6)

missval <- ifelse(!is.na(persons$age) & !is.na(persons$female) & !is.na(persons$smkyears) & !is.na(persons$qtyears) &
                    !is.na(persons$cpd) & !is.na(persons$race) & !is.na(persons$emp) &
                    !is.na(persons$fam.lung.trend) & !is.na(persons$bmi) & !is.na(persons$edu6),0,1)


anlst$lcds1[missval==0] <- lcrisk(persons[missval==0,],1)[,4]/1000
anlst$lcds2[missval==0] <- lcrisk(persons[missval==0,],2)[,4]/1000
anlst$lcds3[missval==0] <- lcrisk(persons[missval==0,],3)[,4]/1000
anlst$lcds4[missval==0] <- lcrisk(persons[missval==0,],4)[,4]/1000
anlst$lcds5[missval==0] <- lcrisk(persons[missval==0,],5)[,4]/1000

z <- subset(anlst,missval==0)
lyg_res <- sapply(1:nrow(z),lyg)
lyg_res <- t(lyg_res)
anlst$ely <- lyg_res[,1]
anlst$ely_ct <- lyg_res[,2]
anlst$lyg <- lyg_res[,3]
save(anlst,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")