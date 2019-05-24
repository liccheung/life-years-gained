rm(list=ls(all=TRUE))
library(survival)
library(survey)

drest <- function(design){
  EO.hat <- svyratio(num=~death_risk, denom=~died, design)
  eo <- rbind(svytable(~died,design)[2],svytotal(~death_risk,design)[1],coef(EO.hat),vcov(EO.hat))
  names(eo) <- c("obs","exp","E/O","Var(E/O)")
  return(eo)
}

foo <- function(x){
  if (is.numeric(nhis$death_risk[[x]])==1){
    return(nhis$death_risk[[x]])
  } else {
    return(as.numeric(NA))
  }
}

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis2006_09_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$death_risk <- as.numeric(sapply(1:nrow(nhis),foo))
nhis$uspstf.eligible <- ifelse(nhis$age>=55 & nhis$age <= 80 & 
                                 (nhis$current==1 | (nhis$former==1&nhis$qtyears<=15)) &
                                 nhis$pkyr.cat>=30,1,0)
nhis$adj.wt <- nhis$wt/4
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80 & !is.na(death_risk))
uspstf <- subset(master, uspstf.eligible==1)
notuspstf <- subset(master, uspstf.eligible==0)
nhis2006 <- subset(master,year==2006)
nhis2007 <- subset(master,year==2007)
nhis2008 <- subset(master,year==2008)
nhis2009 <- subset(master,year==2009)
age1 <- subset(master,age<50)
age2 <- subset(master,age>=50 & age<60)
age3 <- subset(master,age>=60 & age<70)
age4 <- subset(master,age>=70)
male <- subset(master,female==0)
female <- subset(master,female==1)
white <- subset(master,race==0)
black <- subset(master,race==1)
hisp <- subset(master,race==2)
other <- subset(master,race==3)
bmi1 <- subset(master,bmi<=18.5)
bmi2 <- subset(master,bmi>18.5 & bmi<=20)
bmi3 <- subset(master,bmi>20 & bmi<=25)
bmi4 <- subset(master,bmi>25 & bmi<=30)
bmi5 <- subset(master,bmi>30 & bmi<=35)
bmi6 <- subset(master,bmi>35)
current <- subset(master,current==1)
former <- subset(master,former==1)
cur.cpd1 <- subset(master,cpd<20 & current==1)
cur.cpd2 <- subset(master,cpd>=20 & cpd<30 & current==1)
cur.cpd3 <- subset(master,cpd>=30 & cpd<40 & current==1)
cur.cpd4 <- subset(master,cpd>=40 & current==1)
form.cpd1 <- subset(master,cpd<20 & former==1)
form.cpd2 <- subset(master,cpd>=20 & cpd<30 & former==1)
form.cpd3 <- subset(master,cpd>=30 & cpd<40 & former==1)
form.cpd4 <- subset(master,cpd>=40 & former==1)
nocomor <- subset(master,comorbidities==0)
comor1 <- subset(master,comorbidities==1)
comor2 <- subset(master,comorbidities==2)
comor3p <- subset(master,comorbidities>=3)
eo.1 <- cbind(drest(master),drest(uspstf),drest(notuspstf))
eo.yr.1 <- cbind(drest(nhis2006),drest(nhis2007),drest(nhis2008),drest(nhis2009))
eo.age.1 <- cbind(drest(age1),drest(age2),drest(age3),drest(age4))
eo.gender.1 <- cbind(drest(male),drest(female))
eo.race.1 <- cbind(drest(white),drest(black),drest(hisp),drest(other))
eo.bmi.1 <- cbind(drest(bmi1),drest(bmi2),drest(bmi3),drest(bmi4),drest(bmi5),drest(bmi6))
eo.smkstat.1 <- cbind(drest(current),drest(former))
eo.cpd.1 <- cbind(drest(cur.cpd1),drest(cur.cpd2),drest(cur.cpd3),drest(cur.cpd4),
                  drest(form.cpd1),drest(form.cpd2),drest(form.cpd3),drest(form.cpd4))
eo.comor.1 <- cbind(drest(nocomor),drest(comor1),drest(comor2),drest(comor3p))

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis2006_09_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$death_risk <- as.numeric(sapply(1:nrow(nhis),foo))
nhis$uspstf.eligible <- ifelse(nhis$age>=55 & nhis$age <= 80 & 
                                 (nhis$current==1 | (nhis$former==1&nhis$qtyears<=15)) &
                                 nhis$pkyr.cat>=30,1,0)
nhis$adj.wt <- nhis$wt/4
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80 & !is.na(death_risk))
uspstf <- subset(master, uspstf.eligible==1)
notuspstf <- subset(master, uspstf.eligible==0)
nhis2006 <- subset(master,year==2006)
nhis2007 <- subset(master,year==2007)
nhis2008 <- subset(master,year==2008)
nhis2009 <- subset(master,year==2009)
age1 <- subset(master,age<50)
age2 <- subset(master,age>=50 & age<60)
age3 <- subset(master,age>=60 & age<70)
age4 <- subset(master,age>=70)
male <- subset(master,female==0)
female <- subset(master,female==1)
white <- subset(master,race==0)
black <- subset(master,race==1)
hisp <- subset(master,race==2)
other <- subset(master,race==3)
bmi1 <- subset(master,bmi<=18.5)
bmi2 <- subset(master,bmi>18.5 & bmi<=20)
bmi3 <- subset(master,bmi>20 & bmi<=25)
bmi4 <- subset(master,bmi>25 & bmi<=30)
bmi5 <- subset(master,bmi>30 & bmi<=35)
bmi6 <- subset(master,bmi>35)
current <- subset(master,current==1)
former <- subset(master,former==1)
cur.cpd1 <- subset(master,cpd<20 & current==1)
cur.cpd2 <- subset(master,cpd>=20 & cpd<30 & current==1)
cur.cpd3 <- subset(master,cpd>=30 & cpd<40 & current==1)
cur.cpd4 <- subset(master,cpd>=40 & current==1)
form.cpd1 <- subset(master,cpd<20 & former==1)
form.cpd2 <- subset(master,cpd>=20 & cpd<30 & former==1)
form.cpd3 <- subset(master,cpd>=30 & cpd<40 & former==1)
form.cpd4 <- subset(master,cpd>=40 & former==1)
nocomor <- subset(master,comorbidities==0)
comor1 <- subset(master,comorbidities==1)
comor2 <- subset(master,comorbidities==2)
comor3p <- subset(master,comorbidities>=3)
eo.2 <- cbind(drest(master),drest(uspstf),drest(notuspstf))
eo.yr.2 <- cbind(drest(nhis2006),drest(nhis2007),drest(nhis2008),drest(nhis2009))
eo.age.2 <- cbind(drest(age1),drest(age2),drest(age3),drest(age4))
eo.gender.2 <- cbind(drest(male),drest(female))
eo.race.2 <- cbind(drest(white),drest(black),drest(hisp),drest(other))
eo.bmi.2 <- cbind(drest(bmi1),drest(bmi2),drest(bmi3),drest(bmi4),drest(bmi5),drest(bmi6))
eo.smkstat.2 <- cbind(drest(current),drest(former))
eo.cpd.2 <- cbind(drest(cur.cpd1),drest(cur.cpd2),drest(cur.cpd3),drest(cur.cpd4),
                  drest(form.cpd1),drest(form.cpd2),drest(form.cpd3),drest(form.cpd4))
eo.comor.2 <- cbind(drest(nocomor),drest(comor1),drest(comor2),drest(comor3p))

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis2006_09_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$death_risk <- as.numeric(sapply(1:nrow(nhis),foo))
nhis$uspstf.eligible <- ifelse(nhis$age>=55 & nhis$age <= 80 & 
                                 (nhis$current==1 | (nhis$former==1&nhis$qtyears<=15)) &
                                 nhis$pkyr.cat>=30,1,0)
nhis$adj.wt <- nhis$wt/4
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80 & !is.na(death_risk))
uspstf <- subset(master, uspstf.eligible==1)
notuspstf <- subset(master, uspstf.eligible==0)
nhis2006 <- subset(master,year==2006)
nhis2007 <- subset(master,year==2007)
nhis2008 <- subset(master,year==2008)
nhis2009 <- subset(master,year==2009)
age1 <- subset(master,age<50)
age2 <- subset(master,age>=50 & age<60)
age3 <- subset(master,age>=60 & age<70)
age4 <- subset(master,age>=70)
male <- subset(master,female==0)
female <- subset(master,female==1)
white <- subset(master,race==0)
black <- subset(master,race==1)
hisp <- subset(master,race==2)
other <- subset(master,race==3)
bmi1 <- subset(master,bmi<=18.5)
bmi2 <- subset(master,bmi>18.5 & bmi<=20)
bmi3 <- subset(master,bmi>20 & bmi<=25)
bmi4 <- subset(master,bmi>25 & bmi<=30)
bmi5 <- subset(master,bmi>30 & bmi<=35)
bmi6 <- subset(master,bmi>35)
current <- subset(master,current==1)
former <- subset(master,former==1)
cur.cpd1 <- subset(master,cpd<20 & current==1)
cur.cpd2 <- subset(master,cpd>=20 & cpd<30 & current==1)
cur.cpd3 <- subset(master,cpd>=30 & cpd<40 & current==1)
cur.cpd4 <- subset(master,cpd>=40 & current==1)
form.cpd1 <- subset(master,cpd<20 & former==1)
form.cpd2 <- subset(master,cpd>=20 & cpd<30 & former==1)
form.cpd3 <- subset(master,cpd>=30 & cpd<40 & former==1)
form.cpd4 <- subset(master,cpd>=40 & former==1)
nocomor <- subset(master,comorbidities==0)
comor1 <- subset(master,comorbidities==1)
comor2 <- subset(master,comorbidities==2)
comor3p <- subset(master,comorbidities>=3)
eo.3 <- cbind(drest(master),drest(uspstf),drest(notuspstf))
eo.yr.3 <- cbind(drest(nhis2006),drest(nhis2007),drest(nhis2008),drest(nhis2009))
eo.age.3 <- cbind(drest(age1),drest(age2),drest(age3),drest(age4))
eo.gender.3 <- cbind(drest(male),drest(female))
eo.race.3 <- cbind(drest(white),drest(black),drest(hisp),drest(other))
eo.bmi.3 <- cbind(drest(bmi1),drest(bmi2),drest(bmi3),drest(bmi4),drest(bmi5),drest(bmi6))
eo.smkstat.3 <- cbind(drest(current),drest(former))
eo.cpd.3 <- cbind(drest(cur.cpd1),drest(cur.cpd2),drest(cur.cpd3),drest(cur.cpd4),
                  drest(form.cpd1),drest(form.cpd2),drest(form.cpd3),drest(form.cpd4))
eo.comor.3 <- cbind(drest(nocomor),drest(comor1),drest(comor2),drest(comor3p))

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis2006_09_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$death_risk <- as.numeric(sapply(1:nrow(nhis),foo))
nhis$uspstf.eligible <- ifelse(nhis$age>=55 & nhis$age <= 80 & 
                                 (nhis$current==1 | (nhis$former==1&nhis$qtyears<=15)) &
                                 nhis$pkyr.cat>=30,1,0)
nhis$adj.wt <- nhis$wt/4
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80 & !is.na(death_risk))
uspstf <- subset(master, uspstf.eligible==1)
notuspstf <- subset(master, uspstf.eligible==0)
nhis2006 <- subset(master,year==2006)
nhis2007 <- subset(master,year==2007)
nhis2008 <- subset(master,year==2008)
nhis2009 <- subset(master,year==2009)
age1 <- subset(master,age<50)
age2 <- subset(master,age>=50 & age<60)
age3 <- subset(master,age>=60 & age<70)
age4 <- subset(master,age>=70)
male <- subset(master,female==0)
female <- subset(master,female==1)
white <- subset(master,race==0)
black <- subset(master,race==1)
hisp <- subset(master,race==2)
other <- subset(master,race==3)
bmi1 <- subset(master,bmi<=18.5)
bmi2 <- subset(master,bmi>18.5 & bmi<=20)
bmi3 <- subset(master,bmi>20 & bmi<=25)
bmi4 <- subset(master,bmi>25 & bmi<=30)
bmi5 <- subset(master,bmi>30 & bmi<=35)
bmi6 <- subset(master,bmi>35)
current <- subset(master,current==1)
former <- subset(master,former==1)
cur.cpd1 <- subset(master,cpd<20 & current==1)
cur.cpd2 <- subset(master,cpd>=20 & cpd<30 & current==1)
cur.cpd3 <- subset(master,cpd>=30 & cpd<40 & current==1)
cur.cpd4 <- subset(master,cpd>=40 & current==1)
form.cpd1 <- subset(master,cpd<20 & former==1)
form.cpd2 <- subset(master,cpd>=20 & cpd<30 & former==1)
form.cpd3 <- subset(master,cpd>=30 & cpd<40 & former==1)
form.cpd4 <- subset(master,cpd>=40 & former==1)
nocomor <- subset(master,comorbidities==0)
comor1 <- subset(master,comorbidities==1)
comor2 <- subset(master,comorbidities==2)
comor3p <- subset(master,comorbidities>=3)
eo.4 <- cbind(drest(master),drest(uspstf),drest(notuspstf))
eo.yr.4 <- cbind(drest(nhis2006),drest(nhis2007),drest(nhis2008),drest(nhis2009))
eo.age.4 <- cbind(drest(age1),drest(age2),drest(age3),drest(age4))
eo.gender.4 <- cbind(drest(male),drest(female))
eo.race.4 <- cbind(drest(white),drest(black),drest(hisp),drest(other))
eo.bmi.4 <- cbind(drest(bmi1),drest(bmi2),drest(bmi3),drest(bmi4),drest(bmi5),drest(bmi6))
eo.smkstat.4 <- cbind(drest(current),drest(former))
eo.cpd.4 <- cbind(drest(cur.cpd1),drest(cur.cpd2),drest(cur.cpd3),drest(cur.cpd4),
                  drest(form.cpd1),drest(form.cpd2),drest(form.cpd3),drest(form.cpd4))
eo.comor.4 <- cbind(drest(nocomor),drest(comor1),drest(comor2),drest(comor3p))

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis2006_09_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$death_risk <- as.numeric(sapply(1:nrow(nhis),foo))
nhis$uspstf.eligible <- ifelse(nhis$age>=55 & nhis$age <= 80 & 
                                 (nhis$current==1 | (nhis$former==1&nhis$qtyears<=15)) &
                                 nhis$pkyr.cat>=30,1,0)
nhis$adj.wt <- nhis$wt/4
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80 & !is.na(death_risk))
uspstf <- subset(master, uspstf.eligible==1)
notuspstf <- subset(master, uspstf.eligible==0)
nhis2006 <- subset(master,year==2006)
nhis2007 <- subset(master,year==2007)
nhis2008 <- subset(master,year==2008)
nhis2009 <- subset(master,year==2009)
age1 <- subset(master,age<50)
age2 <- subset(master,age>=50 & age<60)
age3 <- subset(master,age>=60 & age<70)
age4 <- subset(master,age>=70)
male <- subset(master,female==0)
female <- subset(master,female==1)
white <- subset(master,race==0)
black <- subset(master,race==1)
hisp <- subset(master,race==2)
other <- subset(master,race==3)
bmi1 <- subset(master,bmi<=18.5)
bmi2 <- subset(master,bmi>18.5 & bmi<=20)
bmi3 <- subset(master,bmi>20 & bmi<=25)
bmi4 <- subset(master,bmi>25 & bmi<=30)
bmi5 <- subset(master,bmi>30 & bmi<=35)
bmi6 <- subset(master,bmi>35)
current <- subset(master,current==1)
former <- subset(master,former==1)
cur.cpd1 <- subset(master,cpd<20 & current==1)
cur.cpd2 <- subset(master,cpd>=20 & cpd<30 & current==1)
cur.cpd3 <- subset(master,cpd>=30 & cpd<40 & current==1)
cur.cpd4 <- subset(master,cpd>=40 & current==1)
form.cpd1 <- subset(master,cpd<20 & former==1)
form.cpd2 <- subset(master,cpd>=20 & cpd<30 & former==1)
form.cpd3 <- subset(master,cpd>=30 & cpd<40 & former==1)
form.cpd4 <- subset(master,cpd>=40 & former==1)
nocomor <- subset(master,comorbidities==0)
comor1 <- subset(master,comorbidities==1)
comor2 <- subset(master,comorbidities==2)
comor3p <- subset(master,comorbidities>=3)
eo.5 <- cbind(drest(master),drest(uspstf),drest(notuspstf))
eo.yr.5 <- cbind(drest(nhis2006),drest(nhis2007),drest(nhis2008),drest(nhis2009))
eo.age.5 <- cbind(drest(age1),drest(age2),drest(age3),drest(age4))
eo.gender.5 <- cbind(drest(male),drest(female))
eo.race.5 <- cbind(drest(white),drest(black),drest(hisp),drest(other))
eo.bmi.5 <- cbind(drest(bmi1),drest(bmi2),drest(bmi3),drest(bmi4),drest(bmi5),drest(bmi6))
eo.smkstat.5 <- cbind(drest(current),drest(former))
eo.cpd.5 <- cbind(drest(cur.cpd1),drest(cur.cpd2),drest(cur.cpd3),drest(cur.cpd4),
                  drest(form.cpd1),drest(form.cpd2),drest(form.cpd3),drest(form.cpd4))
eo.comor.5 <- cbind(drest(nocomor),drest(comor1),drest(comor2),drest(comor3p))

ave.eo <- (eo.1+eo.2+eo.3+eo.4+eo.5)/5
ave.yr.eo <- (eo.yr.1+eo.yr.2+eo.yr.3+eo.yr.4+eo.yr.5)/5  #Better as an internal check
ave.age.eo <- (eo.age.1+eo.age.2+eo.age.3+eo.age.4+eo.age.5)/5
ave.gender.eo <- (eo.gender.1+eo.gender.2+eo.gender.3+eo.gender.4+eo.gender.5)/5
ave.race.eo <- (eo.race.1+eo.race.2+eo.race.3+eo.race.4+eo.race.5)/5
ave.bmi.eo <- (eo.bmi.1+eo.bmi.2+eo.bmi.3+eo.bmi.4+eo.bmi.5)/5
ave.smkstat.eo <- (eo.smkstat.1+eo.smkstat.2+eo.smkstat.3+eo.smkstat.4+eo.smkstat.5)/5
ave.cpd.eo <- (eo.cpd.1+eo.cpd.2+eo.cpd.3+eo.cpd.4+eo.cpd.5)/5
ave.comor.eo <- (eo.comor.1+eo.comor.2+eo.comor.3+eo.comor.4+eo.comor.5)/5

#average over the 5 simulations							
means.bar <- ave.eo[3,]
vars.bar <- rowMeans(vars)

#v.impute is the sample variance in E/O
#V(E/O) = (E/O - E(E/O))^2/n
v.impute <- rowMeans((means[13:18,]-means.bar[13:18])^2/5)

#Calculate confidence intervals for mean E/O
#logv - delta method derivation of log(E(E/O))?
logv <- (vars.bar+(1+1/5)*v.impute)/means.bar[13:18]^2
lower <- exp(log(means.bar[13:18])-1.96*sqrt(logv))
upper <- exp(log(means.bar[13:18])+1.96*sqrt(logv))