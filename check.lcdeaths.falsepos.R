rm(list=ls(all=TRUE))
library("haven")
library(lcrisks)
library(ggplot2)
library(plotrix)
require(gridExtra)
load(file="~/Desktop/Lung cancer/lrisk/other/nlst/anlst_010918.RData")
anlst2 <- subset(anlst,select=c("pid",
                                "comp_death",
                                "comp_respfail",
                                "comp_anaphylaxis",
                                "comp_bronchfist",
                                "comp_cardiac",
                                "comp_stroke",
                                "comp_heartfail",
                                "comp_hemothorax",
                                "comp_myocardial",
                                "comp_respiratory",
                                "comp_wound",
                                "comp_bronchial",
                                "comp_empyema",
                                "comp_injury",
                                "comp_ventilation",
                                "comp_thromboembolic",
                                "comp_chylousfist",
                                "comp_brachial",
                                "comp_collapse",
                                "comp_colon"))
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")
anlst <- merge(anlst,anlst2,by="pid")

anlst$comp <- pmax(anlst$comp_respfail,anlst$comp_anaphylaxis,anlst$comp_bronchfist,anlst$comp_cardiac,
                   anlst$comp_stroke,anlst$comp_heartfail,anlst$comp_hemothorax,anlst$comp_myocardial,
                   anlst$comp_respiratory,anlst$comp_wound,anlst$comp_bronchial,anlst$comp_empyema,
                   anlst$comp_injury,anlst$comp_ventilation,anlst$comp_thromboembolic,anlst$comp_chylousfist,
                   anlst$comp_brachial,anlst$comp_collapse,anlst$comp_colon)
anlst$comorbidities <- anlst$emp+anlst$hypertension+anlst$chd+anlst$angina+anlst$heartattack+
  anlst$heartdisease+anlst$stroke+anlst$diab+anlst$bron+anlst$kidney+
  anlst$liver+anlst$prior.cancer+anlst$speceq

anlst$t0 <- ifelse(is.na(anlst$truefalse_scrnres_ly0)|anlst$truefalse_scrnres_ly0==0,
                   0,anlst$truefalse_scrnres_ly0)  # EXCLUDE NO/INADEQUATE
anlst$t0 <- ifelse(anlst$t0==0,"No Screen",
                   ifelse(anlst$t0<=2,"True Positive",
                          ifelse(anlst$t0==3,"False Positive","Negative")))

anlst$t1 <- ifelse(is.na(anlst$truefalse_scrnres_ly1)|anlst$truefalse_scrnres_ly1==0,
                   0,anlst$truefalse_scrnres_ly1)  # EXCLUDE NO/INADEQUATE
anlst$t1 <- ifelse(anlst$t1==0,"No Screen",
                   ifelse(anlst$t1<=2,"True Positive",
                          ifelse(anlst$t1==3,"False Positive","Negative")))

anlst$t2 <- ifelse(is.na(anlst$truefalse_scrnres_ly2)|anlst$truefalse_scrnres_ly2==0,
                   0,anlst$truefalse_scrnres_ly2)  # EXCLUDE NO/INADEQUATE
anlst$t2 <- ifelse(anlst$t2==0,"No Screen",
                   ifelse(anlst$t2<=2,"True Positive",
                          ifelse(anlst$t2==3,"False Positive","Negative")))

#121 (32.1%) lung cancer deaths among those with t0 positive screens
#64 (0.94%) lung cancer deaths among those with t0 false positives screens
#160 (0.8%) lung cancer deaths among those with t0 negative screens

#85 (30%) lung cancer deaths among those with t1 positive screens
#33 (0.5%) lung cancer deaths among those with t1 false positives screens
#100 (0.56%) lung cancer deaths among those with t1 negative screens

#54 (24.4%) lung cancer deaths among those with t2 positive screens
#17 (0.44%) lung cancer deaths among those with t2 false positive screens
#82 (0.41%) lung cancer deaths among those with t2 negative screens

#Need to look at survival outcome variables to compare across time points
anlst$incidence.years <- anlst$incidence.days/365.25
prev <- subset(anlst,screen_group=="CT" & t0=="True Positive")
fp_incid <- subset(anlst,screen_group=="CT" & t0=="False Positive" & event.lung.cancer==1)
kmm_prev <- survfit(Surv(prev$death.years,prev$lung.cancer.death)~1)
fp_incid$death.years <- fp_incid$death.years - fp_incid$incidence.years
kmm_fp <- survfit(Surv(fp_incid$death.years,fp_incid$lung.cancer.death)~1)
plot(kmm_prev$time,kmm_prev$surv,type="l",lwd=3,xlab="Years",ylab="Survival",main="Survival after lung cancer detection")
lines(kmm_fp$time,kmm_fp$surv,lwd=3,col="blue")
#need survival variables from lung cancer detection