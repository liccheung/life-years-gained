#This mortality model is the best AIC model using NHIS 97-2009 odd years as training data
rm(list=ls())
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_1.RData")

nhis <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt_mort, data=nhis, nest=TRUE)
design <- subset(master, age>=40&age<=80 & analypop==1 & deathage>age)
#All NHIS year transformations have same AIC
#linear cpd and log(pkyr.cat) barely beats out sqrt(cpd) and sqrt(smkyears)/log(smkyears)
morat <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + year + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    emp + hypertension + chd + angina + heartattack + heartdisease + stroke + diab + bron + kidney + liver + speceq + prior.cancer + 
                    I(log(qtyears + 1)) + cpd + I(log(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

morat_nocomorbid <- svycoxph(Surv(age, deathage, died) ~ 
                    female + race + edu6 + year + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    I(log(qtyears + 1)) + cpd + I(log(pkyr.cat)),
                  design=design, data = nhis, model = TRUE, y = TRUE) 

save(morat, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v5.RData")
