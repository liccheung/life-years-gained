rm(list=ls())
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")

nhis <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
design <- subset(master, age>=40)
morat <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + birthyear + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                    I(log(qtyears + 1)) + cpd + I(sqrt(smkyears)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

morat_nocomorbid <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + birthyear + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                      I(log(qtyears + 1)) + cpd + I(sqrt(smkyears)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

save(morat, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v2.RData")
