rm(list=ls(all=TRUE))
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/nlst/anlst_010918.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v6.RData")

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

anlst$lcds0125[missval==0] <- lcrisk(persons[missval==0,],.125)[,4]/1000
anlst$lcds0250[missval==0] <- lcrisk(persons[missval==0,],.250)[,4]/1000
anlst$lcds0375[missval==0] <- lcrisk(persons[missval==0,],.375)[,4]/1000
anlst$lcds0500[missval==0] <- lcrisk(persons[missval==0,],.500)[,4]/1000
anlst$lcds0625[missval==0] <- lcrisk(persons[missval==0,],.625)[,4]/1000
anlst$lcds0750[missval==0] <- lcrisk(persons[missval==0,],.750)[,4]/1000
anlst$lcds0875[missval==0] <- lcrisk(persons[missval==0,],.875)[,4]/1000
anlst$lcds1[missval==0] <- lcrisk(persons[missval==0,],1)[,4]/1000
anlst$lcds1125[missval==0] <- lcrisk(persons[missval==0,],1.125)[,4]/1000
anlst$lcds1250[missval==0] <- lcrisk(persons[missval==0,],1.250)[,4]/1000
anlst$lcds1375[missval==0] <- lcrisk(persons[missval==0,],1.375)[,4]/1000
anlst$lcds1500[missval==0] <- lcrisk(persons[missval==0,],1.500)[,4]/1000
anlst$lcds1625[missval==0] <- lcrisk(persons[missval==0,],1.625)[,4]/1000
anlst$lcds1750[missval==0] <- lcrisk(persons[missval==0,],1.750)[,4]/1000
anlst$lcds1875[missval==0] <- lcrisk(persons[missval==0,],1.875)[,4]/1000
anlst$lcds2[missval==0] <- lcrisk(persons[missval==0,],2)[,4]/1000
anlst$lcds2125[missval==0] <- lcrisk(persons[missval==0,],2.125)[,4]/1000
anlst$lcds2250[missval==0] <- lcrisk(persons[missval==0,],2.250)[,4]/1000
anlst$lcds2375[missval==0] <- lcrisk(persons[missval==0,],2.375)[,4]/1000
anlst$lcds2500[missval==0] <- lcrisk(persons[missval==0,],2.500)[,4]/1000
anlst$lcds2625[missval==0] <- lcrisk(persons[missval==0,],2.625)[,4]/1000
anlst$lcds2750[missval==0] <- lcrisk(persons[missval==0,],2.750)[,4]/1000
anlst$lcds2875[missval==0] <- lcrisk(persons[missval==0,],2.875)[,4]/1000
anlst$lcds3[missval==0] <- lcrisk(persons[missval==0,],3)[,4]/1000
anlst$lcds3125[missval==0] <- lcrisk(persons[missval==0,],3.125)[,4]/1000
anlst$lcds3250[missval==0] <- lcrisk(persons[missval==0,],3.250)[,4]/1000
anlst$lcds3375[missval==0] <- lcrisk(persons[missval==0,],3.375)[,4]/1000
anlst$lcds3500[missval==0] <- lcrisk(persons[missval==0,],3.500)[,4]/1000
anlst$lcds3625[missval==0] <- lcrisk(persons[missval==0,],3.625)[,4]/1000
anlst$lcds3750[missval==0] <- lcrisk(persons[missval==0,],3.750)[,4]/1000
anlst$lcds3875[missval==0] <- lcrisk(persons[missval==0,],3.875)[,4]/1000
anlst$lcds4[missval==0] <- lcrisk(persons[missval==0,],4)[,4]/1000
anlst$lcds4125[missval==0] <- lcrisk(persons[missval==0,],4.125)[,4]/1000
anlst$lcds4250[missval==0] <- lcrisk(persons[missval==0,],4.250)[,4]/1000
anlst$lcds4375[missval==0] <- lcrisk(persons[missval==0,],4.375)[,4]/1000
anlst$lcds4500[missval==0] <- lcrisk(persons[missval==0,],4.500)[,4]/1000
anlst$lcds4625[missval==0] <- lcrisk(persons[missval==0,],4.625)[,4]/1000
anlst$lcds4750[missval==0] <- lcrisk(persons[missval==0,],4.750)[,4]/1000
anlst$lcds4875[missval==0] <- lcrisk(persons[missval==0,],4.875)[,4]/1000
anlst$lcds5[missval==0] <- lcrisk(persons[missval==0,],5)[,4]/1000

z <- subset(anlst,missval==0)
surv <- data.frame()
for (x in 1:nrow(z)){
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  surv <- rbind(surv,c$surv)
}

save(anlst,surv,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod2.RData")