rm(list=ls())
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")

master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
design <- subset(master, age>=40)
morat <- svycoxph(Surv(deathyear, died) ~ 
                    female + race + edu6 + fam.lung.trend + emp + I(bmi <= 18.5) + I(cpd > 20) + 
                    cuts(pkyr.cat, c(0, 30, 40, 50, Inf), simple.labels = FALSE, right = FALSE) +
                    I(log(age)) + I(log(bmi)) + I(log(qtyears + 1)) + I(log(smkyears)), 
                    design=design, data = nhis, model = FALSE, y = FALSE)

lcd <- svycoxph(Surv(deathyear2, lung.cancer.death) ~ 
                    female + race + edu6 + fam.lung.trend + emp + I(bmi <= 18.5) + I(cpd > 20) + 
                    cuts(pkyr.cat, c(0, 30, 40, 50, Inf), simple.labels = FALSE, right = FALSE) +
                    I(log(age)) + I(log(bmi)) + I(log(qtyears + 1)) + I(log(smkyears)), 
                    design=design, data = nhis, model = FALSE, y = FALSE)