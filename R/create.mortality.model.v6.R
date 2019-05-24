#This mortality model is the best AIC model using NHIS 97-2009 odd years as training data
rm(list=ls())
library(survival)
library(survey)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_09_odd_1.RData")

nhis$deathage[nhis$deathage<=nhis$age] <- nhis$age[nhis$deathage<=nhis$age] + .125 #give some follow-up time to those whose deathage is <= current age
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt_mort, data=nhis, nest=TRUE)
design <- subset(master, age>=40&age<=84 & analypop==1 & deathage>age & lung.cancer.before==0)
#All NHIS year transformations have same AIC
#beats sqrt(cpd) by decent margin
morat <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + year + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                    I(log(qtyears + 1)) + I(log(cpd)) + I(sqrt(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 
write.table(summary(morat)$conf.int,"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/HR_v6.csv",sep=",")

morat_nocomorbid <- svycoxph(Surv(age, deathage, died) ~ 
                               female + race + edu6 + year + 
                               I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                               I(log(qtyears + 1)) + I(log(cpd)) + I(sqrt(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 
write.table(summary(morat_nocomorbid)$conf.int,"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/HR_nocomorbid_v6.csv",sep=",")

save(morat, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v6.RData")

#Investigating interactions
morat1 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     I(bmi <= 18.5)*I(log(qtyears + 1)) + 
                     I(bmi > 18.5 & bmi <= 20)*I(log(qtyears + 1)) + 
                     I(bmi > 25 & bmi <= 30)*I(log(qtyears + 1)) + 
                     I(bmi > 30 & bmi <= 35)*I(log(qtyears + 1)) + 
                     I(bmi > 35)*I(log(qtyears + 1)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat2 <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + year + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    I(bmi <= 18.5)*I(sqrt(pkyr.cat)) + 
                    I(bmi > 18.5 & bmi <= 20)*I(sqrt(pkyr.cat)) + 
                    I(bmi > 25 & bmi <= 30)*I(sqrt(pkyr.cat)) + 
                    I(bmi > 30 & bmi <= 35)*I(sqrt(pkyr.cat)) + 
                    I(bmi > 35)*I(sqrt(pkyr.cat)) + 
                    emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                    I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

morat3 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     I(bmi <= 18.5)*I(log(cpd + 1)) + 
                     I(bmi > 18.5 & bmi <= 20)*I(log(cpd + 1)) + 
                     I(bmi > 25 & bmi <= 30)*I(log(cpd + 1)) + 
                     I(bmi > 30 & bmi <= 35)*I(log(cpd + 1)) + 
                     I(bmi > 35)*I(log(cpd + 1)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat4 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     I(bmi <= 18.5)*hypertension + 
                     I(bmi > 18.5 & bmi <= 20)*hypertension + 
                     I(bmi > 25 & bmi <= 30)*hypertension + 
                     I(bmi > 30 & bmi <= 35)*hypertension + 
                     I(bmi > 35)*hypertension + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 


morat5 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     I(bmi <= 18.5)*liver + 
                     I(bmi > 18.5 & bmi <= 20)*liver + 
                     I(bmi > 25 & bmi <= 30)*liver + 
                     I(bmi > 30 & bmi <= 35)*liver + 
                     I(bmi > 35)*liver + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat6 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     I(bmi <= 18.5)*race + 
                     I(bmi > 18.5 & bmi <= 20)*race + 
                     I(bmi > 25 & bmi <= 30)*race + 
                     I(bmi > 30 & bmi <= 35)*race + 
                     I(bmi > 35)*race + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat7 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     female*edu6 + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat8 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     edu6*I(log(cpd + 1)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat9 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     edu6*I(sqrt(pkyr.cat)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat10 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     race*I(log(qtyears + 1)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 

morat11 <- svycoxph(Surv(age, deathage, died) ~ 
                      female + race + edu6 + year + 
                      I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                      year*I(log(cpd + 1))  + 
                      emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                      I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                    design=design, data = nhis, model = TRUE, y = TRUE) 

morat12 <- svycoxph(Surv(age, deathage, died) ~ 
                     female + race + edu6 + year + 
                     I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                     year*I(sqrt(pkyr.cat)) + 
                     emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                     I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                   design=design, data = nhis, model = TRUE, y = TRUE) 
