rm(list=ls())
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")

lyg <- function(x){
  print(x)
  c <- survfit(morat,newdata=nhis[x,],start.time=nhis$age[x])
  ely <- sum(diff(c(nhis$age[x],c$time))*c$surv)  
  c <- survfit(morat,newdata=nhis[x,],start.time=nhis$age[x]+5,stop.time=100)
  lyg <- (1-nhis$lcd1[x])**.8-(1-nhis$lcd1[x]) + 
    (1-nhis$lcd2[x])**.8-(1-nhis$lcd2[x]) +
    (1-nhis$lcd3[x])**.8-(1-nhis$lcd3[x]) +
    (1-nhis$lcd4[x])**.8-(1-nhis$lcd4[x]) +
    ((1-nhis$lcd5[x])**.8-(1-nhis$lcd5[x]))*(1+sum(diff(c(nhis$age[x]+5,c$time))*c$surv))
  res <- c(ely,ely+lyg,lyg)
  return(res)
}

nhis <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
design <- subset(master, age>=40)
morat <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + birthyear + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                    I(log(qtyears + 1)) + cpd + I(sqrt(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

morat_nocomorbid <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + birthyear + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    I(log(qtyears + 1)) + cpd + I(sqrt(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

save(morat, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_1.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_1.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_2.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_2.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_3.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_3.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_4.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_4.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_5.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_5.RData")