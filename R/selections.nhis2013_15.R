#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==0 | (nhis$age>=40 & nhis$age<=80), nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$LCRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$LCRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$C1 <- ifelse(nhis$uspstf.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcrat <- subset(master,lcrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_lcrat <- subset(master,uspstf.eligible&lcrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcrat.eligible&lyg.eligible)
lcrat_lyg <- subset(master,!uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcrat.eligible&!lyg.eligible)
lcrat_only <- subset(master,!uspstf.eligible&lcrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab0.1 <- svytable(~uspstf.eligible|lcrat.eligible|lyg.eligible,master)
tab1.1 <- svytable(~uspstf.eligible+lcrat.eligible,master)
tab2.1 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.1 <- svytable(~lcrat.eligible+lyg.eligible,master)
tab4.1 <- svytable(~uspstf.eligible&lcrat.eligible&lyg.eligible,master)
tab5.1 <- svytable(~uspstf.eligible&lcrat.eligible&!lyg.eligible,master)
tab6.1 <- svytable(~uspstf.eligible&!lcrat.eligible&lyg.eligible,master)
tab7.1 <- svytable(~uspstf.eligible==0&lcrat.eligible&lyg.eligible,master)
tab8.1 <- svytable(~uspstf.eligible&!lcrat.eligible&!lyg.eligible,master)
tab9.1 <- svytable(~uspstf.eligible==0&lcrat.eligible&!lyg.eligible,master)
tab10.1 <- svytable(~uspstf.eligible==0&!lcrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.1 <- svytotal(~LCRAT,master)
lc1.1 <- svytotal(~LCRAT,uspstf)
lc2.1 <- svytotal(~LCRAT,lcrat)
lc3.1 <- svytotal(~LCRAT,lyg)
lc4.1 <- svytotal(~LCRAT,all)
lc5.1 <- svytotal(~LCRAT,uspstf_lcrat)
lc6.1 <- svytotal(~LCRAT,uspstf_lyg)
lc7.1 <- svytotal(~LCRAT,lcrat_lyg)
lc8.1 <- svytotal(~LCRAT,uspstf_only)
lc9.1 <- svytotal(~LCRAT,lcrat_only)
lc10.1 <- svytotal(~LCRAT,lyg_only)

#Number of Lung Cancer Deaths
lcd0.1 <- svytotal(~lcd5,master)
lcd1.1 <- svytotal(~lcd5,uspstf)
lcd2.1 <- svytotal(~lcd5,lcrat)
lcd3.1 <- svytotal(~lcd5,lyg)
lcd4.1 <- svytotal(~lcd5,all)
lcd5.1 <- svytotal(~lcd5,uspstf_lcrat)
lcd6.1 <- svytotal(~lcd5,uspstf_lyg)
lcd7.1 <- svytotal(~lcd5,lcrat_lyg)
lcd8.1 <- svytotal(~lcd5,uspstf_only)
lcd9.1 <- svytotal(~lcd5,lcrat_only)
lcd10.1 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.1 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.1 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.1 <- svytotal(~lcd5,lcrat)*(1-0.796)
lcds3.1 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.1 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.1 <- svytotal(~lcd5,uspstf_lcrat)*(1-0.796)
lcds6.1 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.1 <- svytotal(~lcd5,lcrat_lyg)*(1-0.796)
lcds8.1 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.1 <- svytotal(~lcd5,lcrat_only)*(1-0.796)
lcds10.1 <- svytotal(~lcd5,lyg_only)*(1-0.796)

#Number of Life-years Gained
lyg0.1 <- svytotal(~lyg,master)
lyg1.1 <- svytotal(~lyg,uspstf)
lyg2.1 <- svytotal(~lyg,lcrat)
lyg3.1 <- svytotal(~lyg,lyg)
lyg4.1 <- svytotal(~lyg,all)
lyg5.1 <- svytotal(~lyg,uspstf_lcrat)
lyg6.1 <- svytotal(~lyg,uspstf_lyg)
lyg7.1 <- svytotal(~lyg,lcrat_lyg)
lyg8.1 <- svytotal(~lyg,uspstf_only)
lyg9.1 <- svytotal(~lyg,lcrat_only)
lyg10.1 <- svytotal(~lyg,lyg_only)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCRAT.q.1 <- svyquantile(~LCRAT,any,quantiles)
nhis$LCRATcat <- ifelse(nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[2],1,
                       ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[2] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[3],2,
                              ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[3] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[4],3,
                                     ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[4] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[5],4,5))))

lyg.q.1 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCRAT0.1 <- svymean(~LCRAT,master) 
mLCRAT1.1 <- svymean(~LCRAT,c1)
mLCRAT2.1 <- svymean(~LCRAT,c2)
mLCRAT3.1 <- svymean(~LCRAT,c3)
mlyg0.1 <- svymean(~lyg,master) 
mlyg1.1 <- svymean(~lyg,c1)
mlyg2.1 <- svymean(~lyg,c2)
mlyg3.1 <- svymean(~lyg,c3)

tabLCRAT0.1 <- svytable(~LCRATcat,master)
tabLCRAT1.1 <- svytable(~LCRATcat+C1,master)
tabLCRAT2.1 <- svytable(~LCRATcat+C2,master)
tabLCRAT3.1 <- svytable(~LCRATcat+C3,master)

tablyg0.1 <- svytable(~lygcat,master)
tablyg1.1 <- svytable(~lygcat+C1,master)
tablyg2.1 <- svytable(~lygcat+C2,master)
tablyg3.1 <- svytable(~lygcat+C3,master)

tabg0.1 <- svytable(~female,master)
tabg1.1 <- svytable(~female+C1,master)
tabg2.1 <- svytable(~female+C2,master)
tabg3.1 <- svytable(~female+C3,master)

tabr0.1 <- svytable(~race,master)
tabr1.1 <- svytable(~race+C1,master)
tabr2.1 <- svytable(~race+C2,master)
tabr3.1 <- svytable(~race+C3,master)

meda0.1 <- svyquantile(~age,master,.5)
meda.uspstf.1 <- svyquantile(~age,uspstf,.5)
meda.LCRAT.1 <- svyquantile(~age,lcrat,.5)
meda.lyg.1 <- svyquantile(~age,lyg,.5)
meda1.1 <- svyquantile(~age,c1,.5)
meda2.1 <- svyquantile(~age,c2,.5)
meda3.1 <- svyquantile(~age,c3,.5)

taba0.1 <- svytable(~age.cat,master)
taba1.1 <- svytable(~age.cat+C1,master)
taba2.1 <- svytable(~age.cat+C2,master)
taba3.1 <- svytable(~age.cat+C3,master)

tabc0.1 <- svytable(~comorbidcat,master)
tabc1.1 <- svytable(~comorbidcat+C1,master)
tabc2.1 <- svytable(~comorbidcat+C2,master)
tabc3.1 <- svytable(~comorbidcat+C3,master)

agec2.1 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.1 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==0 | (nhis$age>=40 & nhis$age<=80), nhis$lyg,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$LCRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.2 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$LCRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.2 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$C1 <- ifelse(nhis$uspstf.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcrat <- subset(master,lcrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_lcrat <- subset(master,uspstf.eligible&lcrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcrat.eligible&lyg.eligible)
lcrat_lyg <- subset(master,!uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcrat.eligible&!lyg.eligible)
lcrat_only <- subset(master,!uspstf.eligible&lcrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab0.2 <- svytable(~uspstf.eligible|lcrat.eligible|lyg.eligible,master)
tab1.2 <- svytable(~uspstf.eligible+lcrat.eligible,master)
tab2.2 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.2 <- svytable(~lcrat.eligible+lyg.eligible,master)
tab4.2 <- svytable(~uspstf.eligible&lcrat.eligible&lyg.eligible,master)
tab5.2 <- svytable(~uspstf.eligible&lcrat.eligible&!lyg.eligible,master)
tab6.2 <- svytable(~uspstf.eligible&!lcrat.eligible&lyg.eligible,master)
tab7.2 <- svytable(~uspstf.eligible==0&lcrat.eligible&lyg.eligible,master)
tab8.2 <- svytable(~uspstf.eligible&!lcrat.eligible&!lyg.eligible,master)
tab9.2 <- svytable(~uspstf.eligible==0&lcrat.eligible&!lyg.eligible,master)
tab10.2 <- svytable(~uspstf.eligible==0&!lcrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.2 <- svytotal(~LCRAT,master)
lc1.2 <- svytotal(~LCRAT,uspstf)
lc2.2 <- svytotal(~LCRAT,lcrat)
lc3.2 <- svytotal(~LCRAT,lyg)
lc4.2 <- svytotal(~LCRAT,all)
lc5.2 <- svytotal(~LCRAT,uspstf_lcrat)
lc6.2 <- svytotal(~LCRAT,uspstf_lyg)
lc7.2 <- svytotal(~LCRAT,lcrat_lyg)
lc8.2 <- svytotal(~LCRAT,uspstf_only)
lc9.2 <- svytotal(~LCRAT,lcrat_only)
lc10.2 <- svytotal(~LCRAT,lyg_only)

#Number of Lung Cancer Deaths
lcd0.2 <- svytotal(~lcd5,master)
lcd1.2 <- svytotal(~lcd5,uspstf)
lcd2.2 <- svytotal(~lcd5,lcrat)
lcd3.2 <- svytotal(~lcd5,lyg)
lcd4.2 <- svytotal(~lcd5,all)
lcd5.2 <- svytotal(~lcd5,uspstf_lcrat)
lcd6.2 <- svytotal(~lcd5,uspstf_lyg)
lcd7.2 <- svytotal(~lcd5,lcrat_lyg)
lcd8.2 <- svytotal(~lcd5,uspstf_only)
lcd9.2 <- svytotal(~lcd5,lcrat_only)
lcd10.2 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.2 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.2 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.2 <- svytotal(~lcd5,lcrat)*(1-0.796)
lcds3.2 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.2 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.2 <- svytotal(~lcd5,uspstf_lcrat)*(1-0.796)
lcds6.2 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.2 <- svytotal(~lcd5,lcrat_lyg)*(1-0.796)
lcds8.2 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.2 <- svytotal(~lcd5,lcrat_only)*(1-0.796)
lcds10.2 <- svytotal(~lcd5,lyg_only)*(1-0.796)

#Number of Life-years Gained
lyg0.2 <- svytotal(~lyg,master)
lyg1.2 <- svytotal(~lyg,uspstf)
lyg2.2 <- svytotal(~lyg,lcrat)
lyg3.2 <- svytotal(~lyg,lyg)
lyg4.2 <- svytotal(~lyg,all)
lyg5.2 <- svytotal(~lyg,uspstf_lcrat)
lyg6.2 <- svytotal(~lyg,uspstf_lyg)
lyg7.2 <- svytotal(~lyg,lcrat_lyg)
lyg8.2 <- svytotal(~lyg,uspstf_only)
lyg9.2 <- svytotal(~lyg,lcrat_only)
lyg10.2 <- svytotal(~lyg,lyg_only)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCRAT.q.2 <- svyquantile(~LCRAT,any,quantiles)
nhis$LCRATcat <- ifelse(nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[2],1,
                        ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[2] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[3],2,
                               ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[3] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[4],3,
                                      ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[4] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[5],4,5))))

lyg.q.2 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCRAT0.2 <- svymean(~LCRAT,master) 
mLCRAT1.2 <- svymean(~LCRAT,c1)
mLCRAT2.2 <- svymean(~LCRAT,c2)
mLCRAT3.2 <- svymean(~LCRAT,c3)
mlyg0.2 <- svymean(~lyg,master) 
mlyg1.2 <- svymean(~lyg,c1)
mlyg2.2 <- svymean(~lyg,c2)
mlyg3.2 <- svymean(~lyg,c3)

tabLCRAT0.2 <- svytable(~LCRATcat,master)
tabLCRAT1.2 <- svytable(~LCRATcat+C1,master)
tabLCRAT2.2 <- svytable(~LCRATcat+C2,master)
tabLCRAT3.2 <- svytable(~LCRATcat+C3,master)

tablyg0.2 <- svytable(~lygcat,master)
tablyg1.2 <- svytable(~lygcat+C1,master)
tablyg2.2 <- svytable(~lygcat+C2,master)
tablyg3.2 <- svytable(~lygcat+C3,master)

tabg0.2 <- svytable(~female,master)
tabg1.2 <- svytable(~female+C1,master)
tabg2.2 <- svytable(~female+C2,master)
tabg3.2 <- svytable(~female+C3,master)

tabr0.2 <- svytable(~race,master)
tabr1.2 <- svytable(~race+C1,master)
tabr2.2 <- svytable(~race+C2,master)
tabr3.2 <- svytable(~race+C3,master)

meda0.2 <- svyquantile(~age,master,.5)
meda.uspstf.2 <- svyquantile(~age,uspstf,.5)
meda.LCRAT.2 <- svyquantile(~age,lcrat,.5)
meda.lyg.2 <- svyquantile(~age,lyg,.5)
meda1.2 <- svyquantile(~age,c1,.5)
meda2.2 <- svyquantile(~age,c2,.5)
meda3.2 <- svyquantile(~age,c3,.5)

taba0.2 <- svytable(~age.cat,master)
taba1.2 <- svytable(~age.cat+C1,master)
taba2.2 <- svytable(~age.cat+C2,master)
taba3.2 <- svytable(~age.cat+C3,master)

tabc0.2 <- svytable(~comorbidcat,master)
tabc1.2 <- svytable(~comorbidcat+C1,master)
tabc2.2 <- svytable(~comorbidcat+C2,master)
tabc3.2 <- svytable(~comorbidcat+C3,master)

agec2.2 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.2 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==0 | (nhis$age>=40 & nhis$age<=80), nhis$lyg,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$LCRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.3 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$LCRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.3 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$C1 <- ifelse(nhis$uspstf.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcrat <- subset(master,lcrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_lcrat <- subset(master,uspstf.eligible&lcrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcrat.eligible&lyg.eligible)
lcrat_lyg <- subset(master,!uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcrat.eligible&!lyg.eligible)
lcrat_only <- subset(master,!uspstf.eligible&lcrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab0.3 <- svytable(~uspstf.eligible|lcrat.eligible|lyg.eligible,master)
tab1.3 <- svytable(~uspstf.eligible+lcrat.eligible,master)
tab2.3 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.3 <- svytable(~lcrat.eligible+lyg.eligible,master)
tab4.3 <- svytable(~uspstf.eligible&lcrat.eligible&lyg.eligible,master)
tab5.3 <- svytable(~uspstf.eligible&lcrat.eligible&!lyg.eligible,master)
tab6.3 <- svytable(~uspstf.eligible&!lcrat.eligible&lyg.eligible,master)
tab7.3 <- svytable(~uspstf.eligible==0&lcrat.eligible&lyg.eligible,master)
tab8.3 <- svytable(~uspstf.eligible&!lcrat.eligible&!lyg.eligible,master)
tab9.3 <- svytable(~uspstf.eligible==0&lcrat.eligible&!lyg.eligible,master)
tab10.3 <- svytable(~uspstf.eligible==0&!lcrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.3 <- svytotal(~LCRAT,master)
lc1.3 <- svytotal(~LCRAT,uspstf)
lc2.3 <- svytotal(~LCRAT,lcrat)
lc3.3 <- svytotal(~LCRAT,lyg)
lc4.3 <- svytotal(~LCRAT,all)
lc5.3 <- svytotal(~LCRAT,uspstf_lcrat)
lc6.3 <- svytotal(~LCRAT,uspstf_lyg)
lc7.3 <- svytotal(~LCRAT,lcrat_lyg)
lc8.3 <- svytotal(~LCRAT,uspstf_only)
lc9.3 <- svytotal(~LCRAT,lcrat_only)
lc10.3 <- svytotal(~LCRAT,lyg_only)

#Number of Lung Cancer Deaths
lcd0.3 <- svytotal(~lcd5,master)
lcd1.3 <- svytotal(~lcd5,uspstf)
lcd2.3 <- svytotal(~lcd5,lcrat)
lcd3.3 <- svytotal(~lcd5,lyg)
lcd4.3 <- svytotal(~lcd5,all)
lcd5.3 <- svytotal(~lcd5,uspstf_lcrat)
lcd6.3 <- svytotal(~lcd5,uspstf_lyg)
lcd7.3 <- svytotal(~lcd5,lcrat_lyg)
lcd8.3 <- svytotal(~lcd5,uspstf_only)
lcd9.3 <- svytotal(~lcd5,lcrat_only)
lcd10.3 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.3 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.3 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.3 <- svytotal(~lcd5,lcrat)*(1-0.796)
lcds3.3 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.3 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.3 <- svytotal(~lcd5,uspstf_lcrat)*(1-0.796)
lcds6.3 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.3 <- svytotal(~lcd5,lcrat_lyg)*(1-0.796)
lcds8.3 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.3 <- svytotal(~lcd5,lcrat_only)*(1-0.796)
lcds10.3 <- svytotal(~lcd5,lyg_only)*(1-0.796)

#Number of Life-years Gained
lyg0.3 <- svytotal(~lyg,master)
lyg1.3 <- svytotal(~lyg,uspstf)
lyg2.3 <- svytotal(~lyg,lcrat)
lyg3.3 <- svytotal(~lyg,lyg)
lyg4.3 <- svytotal(~lyg,all)
lyg5.3 <- svytotal(~lyg,uspstf_lcrat)
lyg6.3 <- svytotal(~lyg,uspstf_lyg)
lyg7.3 <- svytotal(~lyg,lcrat_lyg)
lyg8.3 <- svytotal(~lyg,uspstf_only)
lyg9.3 <- svytotal(~lyg,lcrat_only)
lyg10.3 <- svytotal(~lyg,lyg_only)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCRAT.q.3 <- svyquantile(~LCRAT,any,quantiles)
nhis$LCRATcat <- ifelse(nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[2],1,
                        ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[2] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[3],2,
                               ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[3] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[4],3,
                                      ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[4] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[5],4,5))))

lyg.q.3 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCRAT0.3 <- svymean(~LCRAT,master) 
mLCRAT1.3 <- svymean(~LCRAT,c1)
mLCRAT2.3 <- svymean(~LCRAT,c2)
mLCRAT3.3 <- svymean(~LCRAT,c3)
mlyg0.3 <- svymean(~lyg,master) 
mlyg1.3 <- svymean(~lyg,c1)
mlyg2.3 <- svymean(~lyg,c2)
mlyg3.3 <- svymean(~lyg,c3)

tabLCRAT0.3 <- svytable(~LCRATcat,master)
tabLCRAT1.3 <- svytable(~LCRATcat+C1,master)
tabLCRAT2.3 <- svytable(~LCRATcat+C2,master)
tabLCRAT3.3 <- svytable(~LCRATcat+C3,master)

tablyg0.3 <- svytable(~lygcat,master)
tablyg1.3 <- svytable(~lygcat+C1,master)
tablyg2.3 <- svytable(~lygcat+C2,master)
tablyg3.3 <- svytable(~lygcat+C3,master)

tabg0.3 <- svytable(~female,master)
tabg1.3 <- svytable(~female+C1,master)
tabg2.3 <- svytable(~female+C2,master)
tabg3.3 <- svytable(~female+C3,master)

tabr0.3 <- svytable(~race,master)
tabr1.3 <- svytable(~race+C1,master)
tabr2.3 <- svytable(~race+C2,master)
tabr3.3 <- svytable(~race+C3,master)

meda0.3 <- svyquantile(~age,master,.5)
meda.uspstf.3 <- svyquantile(~age,uspstf,.5)
meda.LCRAT.3 <- svyquantile(~age,lcrat,.5)
meda.lyg.3 <- svyquantile(~age,lyg,.5)
meda1.3 <- svyquantile(~age,c1,.5)
meda2.3 <- svyquantile(~age,c2,.5)
meda3.3 <- svyquantile(~age,c3,.5)

taba0.3 <- svytable(~age.cat,master)
taba1.3 <- svytable(~age.cat+C1,master)
taba2.3 <- svytable(~age.cat+C2,master)
taba3.3 <- svytable(~age.cat+C3,master)

tabc0.3 <- svytable(~comorbidcat,master)
tabc1.3 <- svytable(~comorbidcat+C1,master)
tabc2.3 <- svytable(~comorbidcat+C2,master)
tabc3.3 <- svytable(~comorbidcat+C3,master)

agec2.3 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.3 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==0 | (nhis$age>=40 & nhis$age<=80), nhis$lyg,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$LCRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.4 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$LCRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.4 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$C1 <- ifelse(nhis$uspstf.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcrat <- subset(master,lcrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_lcrat <- subset(master,uspstf.eligible&lcrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcrat.eligible&lyg.eligible)
lcrat_lyg <- subset(master,!uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcrat.eligible&!lyg.eligible)
lcrat_only <- subset(master,!uspstf.eligible&lcrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab0.4 <- svytable(~uspstf.eligible|lcrat.eligible|lyg.eligible,master)
tab1.4 <- svytable(~uspstf.eligible+lcrat.eligible,master)
tab2.4 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.4 <- svytable(~lcrat.eligible+lyg.eligible,master)
tab4.4 <- svytable(~uspstf.eligible&lcrat.eligible&lyg.eligible,master)
tab5.4 <- svytable(~uspstf.eligible&lcrat.eligible&!lyg.eligible,master)
tab6.4 <- svytable(~uspstf.eligible&!lcrat.eligible&lyg.eligible,master)
tab7.4 <- svytable(~uspstf.eligible==0&lcrat.eligible&lyg.eligible,master)
tab8.4 <- svytable(~uspstf.eligible&!lcrat.eligible&!lyg.eligible,master)
tab9.4 <- svytable(~uspstf.eligible==0&lcrat.eligible&!lyg.eligible,master)
tab10.4 <- svytable(~uspstf.eligible==0&!lcrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.4 <- svytotal(~LCRAT,master)
lc1.4 <- svytotal(~LCRAT,uspstf)
lc2.4 <- svytotal(~LCRAT,lcrat)
lc3.4 <- svytotal(~LCRAT,lyg)
lc4.4 <- svytotal(~LCRAT,all)
lc5.4 <- svytotal(~LCRAT,uspstf_lcrat)
lc6.4 <- svytotal(~LCRAT,uspstf_lyg)
lc7.4 <- svytotal(~LCRAT,lcrat_lyg)
lc8.4 <- svytotal(~LCRAT,uspstf_only)
lc9.4 <- svytotal(~LCRAT,lcrat_only)
lc10.4 <- svytotal(~LCRAT,lyg_only)

#Number of Lung Cancer Deaths
lcd0.4 <- svytotal(~lcd5,master)
lcd1.4 <- svytotal(~lcd5,uspstf)
lcd2.4 <- svytotal(~lcd5,lcrat)
lcd3.4 <- svytotal(~lcd5,lyg)
lcd4.4 <- svytotal(~lcd5,all)
lcd5.4 <- svytotal(~lcd5,uspstf_lcrat)
lcd6.4 <- svytotal(~lcd5,uspstf_lyg)
lcd7.4 <- svytotal(~lcd5,lcrat_lyg)
lcd8.4 <- svytotal(~lcd5,uspstf_only)
lcd9.4 <- svytotal(~lcd5,lcrat_only)
lcd10.4 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.4 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.4 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.4 <- svytotal(~lcd5,lcrat)*(1-0.796)
lcds3.4 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.4 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.4 <- svytotal(~lcd5,uspstf_lcrat)*(1-0.796)
lcds6.4 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.4 <- svytotal(~lcd5,lcrat_lyg)*(1-0.796)
lcds8.4 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.4 <- svytotal(~lcd5,lcrat_only)*(1-0.796)
lcds10.4 <- svytotal(~lcd5,lyg_only)*(1-0.796)

#Number of Life-years Gained
lyg0.4 <- svytotal(~lyg,master)
lyg1.4 <- svytotal(~lyg,uspstf)
lyg2.4 <- svytotal(~lyg,lcrat)
lyg3.4 <- svytotal(~lyg,lyg)
lyg4.4 <- svytotal(~lyg,all)
lyg5.4 <- svytotal(~lyg,uspstf_lcrat)
lyg6.4 <- svytotal(~lyg,uspstf_lyg)
lyg7.4 <- svytotal(~lyg,lcrat_lyg)
lyg8.4 <- svytotal(~lyg,uspstf_only)
lyg9.4 <- svytotal(~lyg,lcrat_only)
lyg10.4 <- svytotal(~lyg,lyg_only)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCRAT.q.4 <- svyquantile(~LCRAT,any,quantiles)
nhis$LCRATcat <- ifelse(nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[2],1,
                        ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[2] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[3],2,
                               ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[3] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[4],3,
                                      ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[4] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[5],4,5))))

lyg.q.4 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCRAT0.4 <- svymean(~LCRAT,master) 
mLCRAT1.4 <- svymean(~LCRAT,c1)
mLCRAT2.4 <- svymean(~LCRAT,c2)
mLCRAT3.4 <- svymean(~LCRAT,c3)
mlyg0.4 <- svymean(~lyg,master) 
mlyg1.4 <- svymean(~lyg,c1)
mlyg2.4 <- svymean(~lyg,c2)
mlyg3.4 <- svymean(~lyg,c3)

tabLCRAT0.4 <- svytable(~LCRATcat,master)
tabLCRAT1.4 <- svytable(~LCRATcat+C1,master)
tabLCRAT2.4 <- svytable(~LCRATcat+C2,master)
tabLCRAT3.4 <- svytable(~LCRATcat+C3,master)

tablyg0.4 <- svytable(~lygcat,master)
tablyg1.4 <- svytable(~lygcat+C1,master)
tablyg2.4 <- svytable(~lygcat+C2,master)
tablyg3.4 <- svytable(~lygcat+C3,master)

tabg0.4 <- svytable(~female,master)
tabg1.4 <- svytable(~female+C1,master)
tabg2.4 <- svytable(~female+C2,master)
tabg3.4 <- svytable(~female+C3,master)

tabr0.4 <- svytable(~race,master)
tabr1.4 <- svytable(~race+C1,master)
tabr2.4 <- svytable(~race+C2,master)
tabr3.4 <- svytable(~race+C3,master)

meda0.4 <- svyquantile(~age,master,.5)
meda.uspstf.4 <- svyquantile(~age,uspstf,.5)
meda.LCRAT.4 <- svyquantile(~age,lcrat,.5)
meda.lyg.4 <- svyquantile(~age,lyg,.5)
meda1.4 <- svyquantile(~age,c1,.5)
meda2.4 <- svyquantile(~age,c2,.5)
meda3.4 <- svyquantile(~age,c3,.5)

taba0.4 <- svytable(~age.cat,master)
taba1.4 <- svytable(~age.cat+C1,master)
taba2.4 <- svytable(~age.cat+C2,master)
taba3.4 <- svytable(~age.cat+C3,master)

tabc0.4 <- svytable(~comorbidcat,master)
tabc1.4 <- svytable(~comorbidcat+C1,master)
tabc2.4 <- svytable(~comorbidcat+C2,master)
tabc3.4 <- svytable(~comorbidcat+C3,master)

agec2.4 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.4 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==0 | (nhis$age>=40 & nhis$age<=80), nhis$lyg,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$LCRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.5 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$LCRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.5 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$C1 <- ifelse(nhis$uspstf.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcrat <- subset(master,lcrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_lcrat <- subset(master,uspstf.eligible&lcrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcrat.eligible&lyg.eligible)
lcrat_lyg <- subset(master,!uspstf.eligible&lcrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcrat.eligible&!lyg.eligible)
lcrat_only <- subset(master,!uspstf.eligible&lcrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab0.5 <- svytable(~uspstf.eligible|lcrat.eligible|lyg.eligible,master)
tab1.5 <- svytable(~uspstf.eligible+lcrat.eligible,master)
tab2.5 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.5 <- svytable(~lcrat.eligible+lyg.eligible,master)
tab4.5 <- svytable(~uspstf.eligible&lcrat.eligible&lyg.eligible,master)
tab5.5 <- svytable(~uspstf.eligible&lcrat.eligible&!lyg.eligible,master)
tab6.5 <- svytable(~uspstf.eligible&!lcrat.eligible&lyg.eligible,master)
tab7.5 <- svytable(~uspstf.eligible==0&lcrat.eligible&lyg.eligible,master)
tab8.5 <- svytable(~uspstf.eligible&!lcrat.eligible&!lyg.eligible,master)
tab9.5 <- svytable(~uspstf.eligible==0&lcrat.eligible&!lyg.eligible,master)
tab10.5 <- svytable(~uspstf.eligible==0&!lcrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.5 <- svytotal(~LCRAT,master)
lc1.5 <- svytotal(~LCRAT,uspstf)
lc2.5 <- svytotal(~LCRAT,lcrat)
lc3.5 <- svytotal(~LCRAT,lyg)
lc4.5 <- svytotal(~LCRAT,all)
lc5.5 <- svytotal(~LCRAT,uspstf_lcrat)
lc6.5 <- svytotal(~LCRAT,uspstf_lyg)
lc7.5 <- svytotal(~LCRAT,lcrat_lyg)
lc8.5 <- svytotal(~LCRAT,uspstf_only)
lc9.5 <- svytotal(~LCRAT,lcrat_only)
lc10.5 <- svytotal(~LCRAT,lyg_only)

#Number of Lung Cancer Deaths
lcd0.5 <- svytotal(~lcd5,master)
lcd1.5 <- svytotal(~lcd5,uspstf)
lcd2.5 <- svytotal(~lcd5,lcrat)
lcd3.5 <- svytotal(~lcd5,lyg)
lcd4.5 <- svytotal(~lcd5,all)
lcd5.5 <- svytotal(~lcd5,uspstf_lcrat)
lcd6.5 <- svytotal(~lcd5,uspstf_lyg)
lcd7.5 <- svytotal(~lcd5,lcrat_lyg)
lcd8.5 <- svytotal(~lcd5,uspstf_only)
lcd9.5 <- svytotal(~lcd5,lcrat_only)
lcd10.5 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.5 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.5 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.5 <- svytotal(~lcd5,lcrat)*(1-0.796)
lcds3.5 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.5 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.5 <- svytotal(~lcd5,uspstf_lcrat)*(1-0.796)
lcds6.5 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.5 <- svytotal(~lcd5,lcrat_lyg)*(1-0.796)
lcds8.5 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.5 <- svytotal(~lcd5,lcrat_only)*(1-0.796)
lcds10.5 <- svytotal(~lcd5,lyg_only)*(1-0.796)

#Number of Life-years Gained
lyg0.5 <- svytotal(~lyg,master)
lyg1.5 <- svytotal(~lyg,uspstf)
lyg2.5 <- svytotal(~lyg,lcrat)
lyg3.5 <- svytotal(~lyg,lyg)
lyg4.5 <- svytotal(~lyg,all)
lyg5.5 <- svytotal(~lyg,uspstf_lcrat)
lyg6.5 <- svytotal(~lyg,uspstf_lyg)
lyg7.5 <- svytotal(~lyg,lcrat_lyg)
lyg8.5 <- svytotal(~lyg,uspstf_only)
lyg9.5 <- svytotal(~lyg,lcrat_only)
lyg10.5 <- svytotal(~lyg,lyg_only)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCRAT.q.5 <- svyquantile(~LCRAT,any,quantiles)
nhis$LCRATcat <- ifelse(nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[2],1,
                        ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[2] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[3],2,
                               ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[3] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[4],3,
                                      ifelse(nhis$LCRAT>svyquantile(~LCRAT,any,quantiles)[4] & nhis$LCRAT<=svyquantile(~LCRAT,any,quantiles)[5],4,5))))

lyg.q.5 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 80)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCRAT0.5 <- svymean(~LCRAT,master) 
mLCRAT1.5 <- svymean(~LCRAT,c1)
mLCRAT2.5 <- svymean(~LCRAT,c2)
mLCRAT3.5 <- svymean(~LCRAT,c3)
mlyg0.5 <- svymean(~lyg,master) 
mlyg1.5 <- svymean(~lyg,c1)
mlyg2.5 <- svymean(~lyg,c2)
mlyg3.5 <- svymean(~lyg,c3)

tabLCRAT0.5 <- svytable(~LCRATcat,master)
tabLCRAT1.5 <- svytable(~LCRATcat+C1,master)
tabLCRAT2.5 <- svytable(~LCRATcat+C2,master)
tabLCRAT3.5 <- svytable(~LCRATcat+C3,master)

tablyg0.5 <- svytable(~lygcat,master)
tablyg1.5 <- svytable(~lygcat+C1,master)
tablyg2.5 <- svytable(~lygcat+C2,master)
tablyg3.5 <- svytable(~lygcat+C3,master)

tabg0.5 <- svytable(~female,master)
tabg1.5 <- svytable(~female+C1,master)
tabg2.5 <- svytable(~female+C2,master)
tabg3.5 <- svytable(~female+C3,master)

tabr0.5 <- svytable(~race,master)
tabr1.5 <- svytable(~race+C1,master)
tabr2.5 <- svytable(~race+C2,master)
tabr3.5 <- svytable(~race+C3,master)

meda0.5 <- svyquantile(~age,master,.5)
meda.uspstf.5 <- svyquantile(~age,uspstf,.5)
meda.LCRAT.5 <- svyquantile(~age,lcrat,.5)
meda.lyg.5 <- svyquantile(~age,lyg,.5)
meda1.5 <- svyquantile(~age,c1,.5)
meda2.5 <- svyquantile(~age,c2,.5)
meda3.5 <- svyquantile(~age,c3,.5)

taba0.5 <- svytable(~age.cat,master)
taba1.5 <- svytable(~age.cat+C1,master)
taba2.5 <- svytable(~age.cat+C2,master)
taba3.5 <- svytable(~age.cat+C3,master)

tabc0.5 <- svytable(~comorbidcat,master)
tabc1.5 <- svytable(~comorbidcat+C1,master)
tabc2.5 <- svytable(~comorbidcat+C2,master)
tabc3.5 <- svytable(~comorbidcat+C3,master)

agec2.5 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.5 <- svytable(~age+C3,master)[,2]/svytable(~age,master)


#Average values
(LCRAT.cutoff.1+LCRAT.cutoff.2+LCRAT.cutoff.3+LCRAT.cutoff.4+LCRAT.cutoff.5)/5
365.25*(lyg.cutoff.1+lyg.cutoff.2+lyg.cutoff.3+lyg.cutoff.4+lyg.cutoff.5)/5
(tab0.1+tab0.2+tab0.3+tab0.4+tab0.5)/5
(tab1.1+tab1.2+tab1.3+tab1.4+tab1.5)/5
(tab2.1+tab2.2+tab2.3+tab2.4+tab2.5)/5
(tab3.1+tab3.2+tab3.3+tab3.4+tab3.5)/5
(tab4.1+tab4.2+tab4.3+tab4.4+tab4.5)/5
(tab5.1+tab5.2+tab5.3+tab5.4+tab5.5)/5
(tab6.1+tab6.2+tab6.3+tab6.4+tab6.5)/5
(tab7.1+tab7.2+tab7.3+tab7.4+tab7.5)/5
(tab8.1+tab8.2+tab8.3+tab8.4+tab8.5)/5
(tab9.1+tab9.2+tab9.3+tab9.4+tab9.5)/5
(tab10.1+tab10.2+tab10.3+tab10.4+tab10.5)/5
(lc0.1+lc0.2+lc0.3+lc0.4+lc0.5)/(5*643642)
(lc1.1+lc1.2+lc1.3+lc1.4+lc1.5)/(5*643642)
(lc2.1+lc2.2+lc2.3+lc2.4+lc2.5)/(5*643642)
(lc3.1+lc3.2+lc3.3+lc3.4+lc3.5)/(5*643642)
(lc4.1+lc4.2+lc4.3+lc4.4+lc4.5)/(5*643642)
(lc5.1+lc5.2+lc5.3+lc5.4+lc5.5)/(5*643642)
(lc6.1+lc6.2+lc6.3+lc6.4+lc6.5)/(5*643642)
(lc7.1+lc7.2+lc7.3+lc7.4+lc7.5)/(5*643642)
(lc8.1+lc8.2+lc8.3+lc8.4+lc8.5)/(5*643642)
(lc9.1+lc9.2+lc9.3+lc9.4+lc9.5)/(5*643642)
(lc10.1+lc10.2+lc10.3+lc10.4+lc10.5)/(5*643642)
(lc0.1+lc0.2+lc0.3+lc0.4+lc0.5)/5
(lc1.1+lc1.2+lc1.3+lc1.4+lc1.5)/5
(lc2.1+lc2.2+lc2.3+lc2.4+lc2.5)/5
(lc3.1+lc3.2+lc3.3+lc3.4+lc3.5)/5
(lc4.1+lc4.2+lc4.3+lc4.4+lc4.5)/5
(lc5.1+lc5.2+lc5.3+lc5.4+lc5.5)/5
(lc6.1+lc6.2+lc6.3+lc6.4+lc6.5)/5
(lc7.1+lc7.2+lc7.3+lc7.4+lc7.5)/5
(lc8.1+lc8.2+lc8.3+lc8.4+lc8.5)/5
(lc9.1+lc9.2+lc9.3+lc9.4+lc9.5)/5
(lc10.1+lc10.2+lc10.3+lc10.4+lc10.5)/5
(lcd0.1+lcd0.2+lcd0.3+lcd0.4+lcd0.5)/5
(lcd1.1+lcd1.2+lcd1.3+lcd1.4+lcd1.5)/5
(lcd2.1+lcd2.2+lcd2.3+lcd2.4+lcd2.5)/5
(lcd3.1+lcd3.2+lcd3.3+lcd3.4+lcd3.5)/5
(lcd4.1+lcd4.2+lcd4.3+lcd4.4+lcd4.5)/5
(lcd5.1+lcd5.2+lcd5.3+lcd5.4+lcd5.5)/5
(lcd6.1+lcd6.2+lcd6.3+lcd6.4+lcd6.5)/5
(lcd7.1+lcd7.2+lcd7.3+lcd7.4+lcd7.5)/5
(lcd8.1+lcd8.2+lcd8.3+lcd8.4+lcd8.5)/5
(lcd9.1+lcd9.2+lcd9.3+lcd9.4+lcd9.5)/5
(lcd10.1+lcd10.2+lcd10.3+lcd10.4+lcd10.5)/5
(lcds0.1+lcds0.2+lcds0.3+lcds0.4+lcds0.5)/5
(lcds1.1+lcds1.2+lcds1.3+lcds1.4+lcds1.5)/5
(lcds2.1+lcds2.2+lcds2.3+lcds2.4+lcds2.5)/5
(lcds3.1+lcds3.2+lcds3.3+lcds3.4+lcds3.5)/5
(lcds4.1+lcds4.2+lcds4.3+lcds4.4+lcds4.5)/5
(lcds5.1+lcds5.2+lcds5.3+lcds5.4+lcds5.5)/5
(lcds6.1+lcds6.2+lcds6.3+lcds6.4+lcds6.5)/5
(lcds7.1+lcds7.2+lcds7.3+lcds7.4+lcds7.5)/5
(lcds8.1+lcds8.2+lcds8.3+lcds8.4+lcds8.5)/5
(lcds9.1+lcds9.2+lcds9.3+lcds9.4+lcds9.5)/5
(lcds10.1+lcds10.2+lcds10.3+lcds10.4+lcds10.5)/5
(lyg0.1+lyg0.2+lyg0.3+lyg0.4+lyg0.5)/5
(lyg1.1+lyg1.2+lyg1.3+lyg1.4+lyg1.5)/5
(lyg2.1+lyg2.2+lyg2.3+lyg2.4+lyg2.5)/5
(lyg3.1+lyg3.2+lyg3.3+lyg3.4+lyg3.5)/5
(lyg4.1+lyg4.2+lyg4.3+lyg4.4+lyg4.5)/5
(lyg5.1+lyg5.2+lyg5.3+lyg5.4+lyg5.5)/5
(lyg6.1+lyg6.2+lyg6.3+lyg6.4+lyg6.5)/5
(lyg7.1+lyg7.2+lyg7.3+lyg7.4+lyg7.5)/5
(lyg8.1+lyg8.2+lyg8.3+lyg8.4+lyg8.5)/5
(lyg9.1+lyg9.2+lyg9.3+lyg9.4+lyg9.5)/5
(lyg10.1+lyg10.2+lyg10.3+lyg10.4+lyg10.5)/5
agec2 <- (agec2.1+agec2.2+agec2.3+agec2.4+agec2.5)/5
agec3 <- (agec3.1+agec3.2+agec3.3+agec3.4+agec3.5)/5


#overall
100*(LCRAT.q.1+LCRAT.q.2+LCRAT.q.3+LCRAT.q.4+LCRAT.q.5)/5
365.25*(lyg.q.1+lyg.q.2+lyg.q.3+lyg.q.4+lyg.q.5)/5
(mLCRAT0.1+mLCRAT0.2+mLCRAT0.3+mLCRAT0.4+mLCRAT0.5)/5
365.25*(mlyg0.1+mlyg0.2+mlyg0.3+mlyg0.4+mlyg0.5)/5
(tabLCRAT0.1+tabLCRAT0.2+tabLCRAT0.3+tabLCRAT0.4+tabLCRAT0.5)/5
(tablyg0.1+tablyg0.2+tablyg0.3+tablyg0.4+tablyg0.5)/5
(tabg0.1+tabg0.2+tabg0.3+tabg0.4+tabg0.5)/5
(tabr0.1+tabr0.2+tabr0.3+tabr0.4+tabr0.5)/5
(meda0.1+meda0.2+meda0.3+meda0.4+meda0.5)/5 #median age
(meda.uspstf.1+meda.uspstf.2+meda.uspstf.3+meda.uspstf.4+meda.uspstf.5)/5
(meda.LCRAT.1+meda.LCRAT.2+meda.LCRAT.3+meda.LCRAT.4+meda.LCRAT.5)/5
(meda.lyg.1+meda.lyg.2+meda.lyg.3+meda.lyg.4+meda.lyg.5)/5
(meda1.1+meda1.2+meda1.3+meda1.4+meda1.5)/5
(meda2.1+meda2.2+meda2.3+meda2.4+meda2.5)/5
(meda3.1+meda3.2+meda3.3+meda3.4+meda3.5)/5
(taba0.1+taba0.2+taba0.3+taba0.4+taba0.5)/5 #age categories
(tabc0.1+tabc0.2+tabc0.3+tabc0.4+tabc0.5)/5


#Revised to match new 6 columns
tab4 <- rbind(100*c((mLCRAT1.1+mLCRAT1.2+mLCRAT1.3+mLCRAT1.4+mLCRAT1.5)/5,
(mLCRAT2.1+mLCRAT2.2+mLCRAT2.3+mLCRAT2.4+mLCRAT2.5)/5,
(mLCRAT3.1+mLCRAT3.2+mLCRAT3.3+mLCRAT3.4+mLCRAT3.5)/5),
c(365.25*(mlyg1.1+mlyg1.2+mlyg1.3+mlyg1.4+mlyg1.5)/5,
365.25*(mlyg2.1+mlyg2.2+mlyg2.3+mlyg2.4+mlyg2.5)/5,
365.25*(mlyg3.1+mlyg3.2+mlyg3.3+mlyg3.4+mlyg3.5)/5),
cbind(((tabLCRAT1.1+tabLCRAT1.2+tabLCRAT1.3+tabLCRAT1.4+tabLCRAT1.5)/5)[,2], #LCRAT5 overall
((tabLCRAT2.1+tabLCRAT2.2+tabLCRAT2.3+tabLCRAT2.4+tabLCRAT2.5)/5)[,2],
((tabLCRAT3.1+tabLCRAT3.2+tabLCRAT3.3+tabLCRAT3.4+tabLCRAT3.5)/5)[,2])/1000,
cbind(((tablyg1.1+tablyg1.2+tablyg1.3+tablyg1.4+tablyg1.5)/5)[,2], #lyg overall
((tablyg2.1+tablyg2.2+tablyg2.3+tablyg2.4+tablyg2.5)/5)[,2],
((tablyg3.1+tablyg3.2+tablyg3.3+tablyg3.4+tablyg3.5)/5)[,2])/1000,
cbind(((tabg1.1+tabg1.2+tabg1.3+tabg1.4+tabg1.5)/5)[,2],  #gender
((tabg2.1+tabg2.2+tabg2.3+tabg2.4+tabg2.5)/5)[,2],
((tabg3.1+tabg3.2+tabg3.3+tabg3.4+tabg3.5)/5)[,2])/1000,
cbind(((tabr1.1+tabr1.2+tabr1.3+tabr1.4+tabr1.5)/5)[,2],  #race
((tabr2.1+tabr2.2+tabr2.3+tabr2.4+tabr2.5)/5)[,2],
((tabr3.1+tabr3.2+tabr3.3+tabr3.4+tabr3.5)/5)[,2])/1000,
cbind(((taba1.1+taba1.2+taba1.3+taba1.4+taba1.5)/5)[,2],
((taba2.1+taba2.2+taba2.3+taba2.4+taba2.5)/5)[,2],
((taba3.1+taba3.2+taba3.3+taba3.4+taba3.5)/5)[,2])/1000,
cbind(((tabc1.1+tabc1.2+tabc1.3+tabc1.4+tabc1.5)/5)[,2], #comorbidities
((tabc2.1+tabc2.2+tabc2.3+tabc2.4+tabc2.5)/5)[,2],
((tabc3.1+tabc3.2+tabc3.3+tabc3.4+tabc3.5)/5)[,2])/1000)

colnames(tab4) <- c("USPSTF", "LCRAT", "LYG")
rownames(tab4) <- c("Mean LCRAT","Mean LDG", 
                    "Q1 LCRAT", "Q2 LCRAT", "Q3 LCRAT", "Q4 LCRAT", "Q5 LCRAT",
                    "Q1 LDG", "Q2 LDG", "Q3 LDG", "Q4 LDG", "Q5 LDG",
                    "Males","Females",
                    "Whites","Blacks","Hispanics","Other",
                    "40-49","50-59","60-69","70-80",
                    "No comorbidities","1 comorbidities","2 comordibities","3+ comordibities")

write.table(round(tab4,2),"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/screen_selections_margins.csv",sep=",")