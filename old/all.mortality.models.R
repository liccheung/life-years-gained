# This program was used to find the best AIC model for overall mortality
# Models were fit to NHIS 1997-2001 data
# Expected-observed based on NHIS 2002-2009

rm(list=ls(all=TRUE))

# Set up
library(survival)
library(survey)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")

nhis <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
design <- subset(master, age>=40)

# PART 1: LUNG CANCER DEATH RISK MODELS
modcall <- function(formula) {                         
   cox.mort <- svycoxph(formula, design = design, data = nhis, model = TRUE, y= TRUE)
   sum.mod <- data.frame(as.character(formula)[3],cox.mort$ll[2],2*length(coef(cox.mort)) - 2*cox.mort$ll[2],cox.mort$concordance[1]/(cox.mort$concordance[1]+cox.mort$concordance[2]),row.names=NULL,check.names=FALSE)
   return(sum.mod)
}

modsum <- data.frame()

#USE SAS TO WRITE THE CALLS
modsum <- rbind(modsum,modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+smkyears))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+smkyears))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+smkyears))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+smkyears))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+smkyears))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+smkyears))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+smkyears))))


modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(log(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(log(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(smkyears^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(smkyears^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(smkyears^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(smkyears^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(smkyears^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(smkyears^2)))))



modsum <- t(cbind(t(modsum),modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(smkyears))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(smkyears))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(smkyears))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(smkyears))))))






modsum <- t(cbind(t(modsum),modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+pkyr.cat))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+pkyr.cat))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+pkyr.cat))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+pkyr.cat))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+pkyr.cat))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+pkyr.cat))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(log(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(log(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(log(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(pkyr.cat^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(pkyr.cat^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(pkyr.cat^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(pkyr.cat^2)))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(pkyr.cat^2)))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(pkyr.cat^2)))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+qtyears+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+bmi+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+qtyears+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(log(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+qtyears+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(bmi^2)+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))



modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+qtyears+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(log(qtyears+1))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(qtyears^2)+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))

modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+cpd+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(log(cpd+1))+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(cpd^2)+I(sqrt(pkyr.cat))))))
modsum <- t(cbind(t(modsum),t(modcall(Surv(age, deathage, died) ~ female+race+edu6+I(bmi <= 18.5)+birthyear+emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer+I(sqrt(bmi))+I(sqrt(qtyears))+I(sqrt(cpd))+I(sqrt(pkyr.cat))))))


modsum <- modsum[order(as.numeric(modsum[,3])),]


rownames(modsum) <- c()
colnames(modsum) <- c("Model","Log-likelihood","AIC","c-stat")
filename <- "~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/all.mortality.csv"
write.table(modsum,filename, sep=",", row.names=FALSE)


