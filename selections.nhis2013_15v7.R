#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

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
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_lcdrat <- subset(master,uspstf.eligible&lcdrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcdrat.eligible&lyg.eligible)
lcdrat_lyg <- subset(master,!uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcdrat.eligible&!lyg.eligible)
lcdrat_only <- subset(master,!uspstf.eligible&lcdrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcdrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)

#Output
#Differences in Selection
tab0.1 <- svytable(~uspstf.eligible|lcdrat.eligible|lyg.eligible,master)
tab1.1 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.1 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.1 <- svytable(~lcdrat.eligible+lyg.eligible,master)
tab4.1 <- svytable(~uspstf.eligible&lcdrat.eligible&lyg.eligible,master)
tab5.1 <- svytable(~uspstf.eligible&lcdrat.eligible&!lyg.eligible,master)
tab6.1 <- svytable(~uspstf.eligible&!lcdrat.eligible&lyg.eligible,master)
tab7.1 <- svytable(~uspstf.eligible==0&lcdrat.eligible&lyg.eligible,master)
tab8.1 <- svytable(~uspstf.eligible&!lcdrat.eligible&!lyg.eligible,master)
tab9.1 <- svytable(~uspstf.eligible==0&lcdrat.eligible&!lyg.eligible,master)
tab10.1 <- svytable(~uspstf.eligible==0&!lcdrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.1 <- svytotal(~LCRAT,master)
lc1.1 <- svytotal(~LCRAT,uspstf)
lc2.1 <- svytotal(~LCRAT,lcdrat)
lc3.1 <- svytotal(~LCRAT,lyg)
lc4.1 <- svytotal(~LCRAT,all)
lc5.1 <- svytotal(~LCRAT,uspstf_lcdrat)
lc6.1 <- svytotal(~LCRAT,uspstf_lyg)
lc7.1 <- svytotal(~LCRAT,lcdrat_lyg)
lc8.1 <- svytotal(~LCRAT,uspstf_only)
lc9.1 <- svytotal(~LCRAT,lcdrat_only)
lc10.1 <- svytotal(~LCRAT,lyg_only)
lc11.1 <- svytotal(~LCRAT,lcdrat_notlyg)
lc12.1 <- svytotal(~LCRAT,lyg_notlcdrat)

#Number of Lung Cancer Deaths
lcd0.1 <- svytotal(~lcd5,master)
lcd1.1 <- svytotal(~lcd5,uspstf)
lcd2.1 <- svytotal(~lcd5,lcdrat)
lcd3.1 <- svytotal(~lcd5,lyg)
lcd4.1 <- svytotal(~lcd5,all)
lcd5.1 <- svytotal(~lcd5,uspstf_lcdrat)
lcd6.1 <- svytotal(~lcd5,uspstf_lyg)
lcd7.1 <- svytotal(~lcd5,lcdrat_lyg)
lcd8.1 <- svytotal(~lcd5,uspstf_only)
lcd9.1 <- svytotal(~lcd5,lcdrat_only)
lcd10.1 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.1 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.1 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.1 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.1 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.1 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.1 <- svytotal(~lcd5,uspstf_lcdrat)*(1-0.796)
lcds6.1 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.1 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds8.1 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.1 <- svytotal(~lcd5,lcdrat_only)*(1-0.796)
lcds10.1 <- svytotal(~lcd5,lyg_only)*(1-0.796)
lcds11.1 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds12.1 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)

#Number of Life-years Gained
lyg0.1 <- svytotal(~lyg,master)
lyg1.1 <- svytotal(~lyg,uspstf)
lyg2.1 <- svytotal(~lyg,lcdrat)
lyg3.1 <- svytotal(~lyg,lyg)
lyg4.1 <- svytotal(~lyg,all)
lyg5.1 <- svytotal(~lyg,uspstf_lcdrat)
lyg6.1 <- svytotal(~lyg,uspstf_lyg)
lyg7.1 <- svytotal(~lyg,lcdrat_lyg)
lyg8.1 <- svytotal(~lyg,uspstf_only)
lyg9.1 <- svytotal(~lyg,lcdrat_only)
lyg10.1 <- svytotal(~lyg,lyg_only)
lyg11.1 <- svytotal(~lyg,lcdrat_notlyg)
lyg12.1 <- svytotal(~lyg,lyg_notlcdrat)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCDRAT.q.1 <- svyquantile(~lcd5,any,quantiles)
nhis$LCDRATcat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
                       ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[2] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[3],2,
                              ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[3] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[4],3,
                                     ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[4] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[5],4,5))))

lyg.q.1 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcdrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCDRAT0.1 <- svymean(~lcd5,master) 
mLCDRAT1.1 <- svymean(~lcd5,c1)
mLCDRAT2.1 <- svymean(~lcd5,c2)
mLCDRAT3.1 <- svymean(~lcd5,c3)
mlyg0.1 <- svymean(~lyg,master) 
mlyg1.1 <- svymean(~lyg,c1)
mlyg2.1 <- svymean(~lyg,c2)
mlyg3.1 <- svymean(~lyg,c3)

tabLCDRAT0.1 <- svytable(~LCDRATcat,master)
tabLCDRAT1.1 <- svytable(~LCDRATcat+C1,master)
tabLCDRAT2.1 <- svytable(~LCDRATcat+C2,master)
tabLCDRAT3.1 <- svytable(~LCDRATcat+C3,master)

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
meda.LCDRAT.1 <- svyquantile(~age,lcdrat,.5)
meda.lyg.1 <- svyquantile(~age,lyg,.5)
meana1.1 <- svymean(~age,c1)
meana2.1 <- svymean(~age,c2)
meana3.1 <- svymean(~age,c3)

taba0.1 <- svytable(~age.cat,master)
taba1.1 <- svytable(~age.cat+C1,master)
taba2.1 <- svytable(~age.cat+C2,master)
taba3.1 <- svytable(~age.cat+C3,master)

tabs0.1 <- svytable(~current,master)
tabs1.1 <- svytable(~current+C1,master)
tabs2.1 <- svytable(~current+C2,master)
tabs3.1 <- svytable(~current+C3,master)

meanc0.1 <- svymean(~comorbidities,master)
meanc1.1 <- svymean(~comorbidities,c1)
meanc2.1 <- svymean(~comorbidities,c2)
meanc3.1 <- svymean(~comorbidities,c3)

tabc0.1 <- svytable(~comorbidcat,master)
tabc1.1 <- svytable(~comorbidcat+C1,master)
tabc2.1 <- svytable(~comorbidcat+C2,master)
tabc3.1 <- svytable(~comorbidcat+C3,master)

tabemp0.1 <- svytable(~emp,master)
tabemp1.1 <- svytable(~emp+C1,master)
tabemp2.1 <- svytable(~emp+C2,master)
tabemp3.1 <- svytable(~emp+C3,master)

tabhype0.1 <- svytable(~hypertension,master)
tabhype1.1 <- svytable(~hypertension+C1,master)
tabhype2.1 <- svytable(~hypertension+C2,master)
tabhype3.1 <- svytable(~hypertension+C3,master)

tabchd0.1 <- svytable(~chd,master)
tabchd1.1 <- svytable(~chd+C1,master)
tabchd2.1 <- svytable(~chd+C2,master)
tabchd3.1 <- svytable(~chd+C3,master)

tabang0.1 <- svytable(~angina,master)
tabang1.1 <- svytable(~angina+C1,master)
tabang2.1 <- svytable(~angina+C2,master)
tabang3.1 <- svytable(~angina+C3,master)

tabha0.1 <- svytable(~heartattack,master)
tabha1.1 <- svytable(~heartattack+C1,master)
tabha2.1 <- svytable(~heartattack+C2,master)
tabha3.1 <- svytable(~heartattack+C3,master)

tabhd0.1 <- svytable(~heartdisease,master)
tabhd1.1 <- svytable(~heartdisease+C1,master)
tabhd2.1 <- svytable(~heartdisease+C2,master)
tabhd3.1 <- svytable(~heartdisease+C3,master)

tabstroke0.1 <- svytable(~stroke,master)
tabstroke1.1 <- svytable(~stroke+C1,master)
tabstroke2.1 <- svytable(~stroke+C2,master)
tabstroke3.1 <- svytable(~stroke+C3,master)

tabdiab0.1 <- svytable(~diab,master)
tabdiab1.1 <- svytable(~diab+C1,master)
tabdiab2.1 <- svytable(~diab+C2,master)
tabdiab3.1 <- svytable(~diab+C3,master)

tabbron0.1 <- svytable(~bron,master)
tabbron1.1 <- svytable(~bron+C1,master)
tabbron2.1 <- svytable(~bron+C2,master)
tabbron3.1 <- svytable(~bron+C3,master)

tabkid0.1 <- svytable(~kidney,master)
tabkid1.1 <- svytable(~kidney+C1,master)
tabkid2.1 <- svytable(~kidney+C2,master)
tabkid3.1 <- svytable(~kidney+C3,master)

tabliv0.1 <- svytable(~liver,master)
tabliv1.1 <- svytable(~liver+C1,master)
tabliv2.1 <- svytable(~liver+C2,master)
tabliv3.1 <- svytable(~liver+C3,master)

tabpc0.1 <- svytable(~prior.cancer,master)
tabpc1.1 <- svytable(~prior.cancer+C1,master)
tabpc2.1 <- svytable(~prior.cancer+C2,master)
tabpc3.1 <- svytable(~prior.cancer+C3,master)

tabse0.1 <- svytable(~speceq,master)
tabse1.1 <- svytable(~speceq+C1,master)
tabse2.1 <- svytable(~speceq+C2,master)
tabse3.1 <- svytable(~speceq+C3,master)

agec2.1 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.1 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

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
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.2 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_lcdrat <- subset(master,uspstf.eligible&lcdrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcdrat.eligible&lyg.eligible)
lcdrat_lyg <- subset(master,!uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcdrat.eligible&!lyg.eligible)
lcdrat_only <- subset(master,!uspstf.eligible&lcdrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcdrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)

#Output
#Differences in Selection
tab0.2 <- svytable(~uspstf.eligible|lcdrat.eligible|lyg.eligible,master)
tab1.2 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.2 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.2 <- svytable(~lcdrat.eligible+lyg.eligible,master)
tab4.2 <- svytable(~uspstf.eligible&lcdrat.eligible&lyg.eligible,master)
tab5.2 <- svytable(~uspstf.eligible&lcdrat.eligible&!lyg.eligible,master)
tab6.2 <- svytable(~uspstf.eligible&!lcdrat.eligible&lyg.eligible,master)
tab7.2 <- svytable(~uspstf.eligible==0&lcdrat.eligible&lyg.eligible,master)
tab8.2 <- svytable(~uspstf.eligible&!lcdrat.eligible&!lyg.eligible,master)
tab9.2 <- svytable(~uspstf.eligible==0&lcdrat.eligible&!lyg.eligible,master)
tab10.2 <- svytable(~uspstf.eligible==0&!lcdrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.2 <- svytotal(~LCRAT,master)
lc1.2 <- svytotal(~LCRAT,uspstf)
lc2.2 <- svytotal(~LCRAT,lcdrat)
lc3.2 <- svytotal(~LCRAT,lyg)
lc4.2 <- svytotal(~LCRAT,all)
lc5.2 <- svytotal(~LCRAT,uspstf_lcdrat)
lc6.2 <- svytotal(~LCRAT,uspstf_lyg)
lc7.2 <- svytotal(~LCRAT,lcdrat_lyg)
lc8.2 <- svytotal(~LCRAT,uspstf_only)
lc9.2 <- svytotal(~LCRAT,lcdrat_only)
lc10.2 <- svytotal(~LCRAT,lyg_only)
lc11.2 <- svytotal(~LCRAT,lcdrat_notlyg)
lc12.2 <- svytotal(~LCRAT,lyg_notlcdrat)

#Number of Lung Cancer Deaths
lcd0.2 <- svytotal(~lcd5,master)
lcd1.2 <- svytotal(~lcd5,uspstf)
lcd2.2 <- svytotal(~lcd5,lcdrat)
lcd3.2 <- svytotal(~lcd5,lyg)
lcd4.2 <- svytotal(~lcd5,all)
lcd5.2 <- svytotal(~lcd5,uspstf_lcdrat)
lcd6.2 <- svytotal(~lcd5,uspstf_lyg)
lcd7.2 <- svytotal(~lcd5,lcdrat_lyg)
lcd8.2 <- svytotal(~lcd5,uspstf_only)
lcd9.2 <- svytotal(~lcd5,lcdrat_only)
lcd10.2 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.2 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.2 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.2 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.2 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.2 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.2 <- svytotal(~lcd5,uspstf_lcdrat)*(1-0.796)
lcds6.2 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.2 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds8.2 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.2 <- svytotal(~lcd5,lcdrat_only)*(1-0.796)
lcds10.2 <- svytotal(~lcd5,lyg_only)*(1-0.796)
lcds11.2 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds12.2 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)

#Number of Life-years Gained
lyg0.2 <- svytotal(~lyg,master)
lyg1.2 <- svytotal(~lyg,uspstf)
lyg2.2 <- svytotal(~lyg,lcdrat)
lyg3.2 <- svytotal(~lyg,lyg)
lyg4.2 <- svytotal(~lyg,all)
lyg5.2 <- svytotal(~lyg,uspstf_lcdrat)
lyg6.2 <- svytotal(~lyg,uspstf_lyg)
lyg7.2 <- svytotal(~lyg,lcdrat_lyg)
lyg8.2 <- svytotal(~lyg,uspstf_only)
lyg9.2 <- svytotal(~lyg,lcdrat_only)
lyg10.2 <- svytotal(~lyg,lyg_only)
lyg11.2 <- svytotal(~lyg,lcdrat_notlyg)
lyg12.2 <- svytotal(~lyg,lyg_notlcdrat)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCDRAT.q.2 <- svyquantile(~lcd5,any,quantiles)
nhis$LCDRATcat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
                        ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[2] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[3],2,
                               ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[3] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[4],3,
                                      ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[4] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[5],4,5))))

lyg.q.2 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcdrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCDRAT0.2 <- svymean(~lcd5,master) 
mLCDRAT1.2 <- svymean(~lcd5,c1)
mLCDRAT2.2 <- svymean(~lcd5,c2)
mLCDRAT3.2 <- svymean(~lcd5,c3)
mlyg0.2 <- svymean(~lyg,master) 
mlyg1.2 <- svymean(~lyg,c1)
mlyg2.2 <- svymean(~lyg,c2)
mlyg3.2 <- svymean(~lyg,c3)

tabLCDRAT0.2 <- svytable(~LCDRATcat,master)
tabLCDRAT1.2 <- svytable(~LCDRATcat+C1,master)
tabLCDRAT2.2 <- svytable(~LCDRATcat+C2,master)
tabLCDRAT3.2 <- svytable(~LCDRATcat+C3,master)

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
meda.LCDRAT.2 <- svyquantile(~age,lcdrat,.5)
meda.lyg.2 <- svyquantile(~age,lyg,.5)
meana1.2 <- svymean(~age,c1)
meana2.2 <- svymean(~age,c2)
meana3.2 <- svymean(~age,c3)

taba0.2 <- svytable(~age.cat,master)
taba1.2 <- svytable(~age.cat+C1,master)
taba2.2 <- svytable(~age.cat+C2,master)
taba3.2 <- svytable(~age.cat+C3,master)

tabs0.2 <- svytable(~current,master)
tabs1.2 <- svytable(~current+C1,master)
tabs2.2 <- svytable(~current+C2,master)
tabs3.2 <- svytable(~current+C3,master)

meanc0.2 <- svymean(~comorbidities,master)
meanc1.2 <- svymean(~comorbidities,c1)
meanc2.2 <- svymean(~comorbidities,c2)
meanc3.2 <- svymean(~comorbidities,c3)

tabc0.2 <- svytable(~comorbidcat,master)
tabc1.2 <- svytable(~comorbidcat+C1,master)
tabc2.2 <- svytable(~comorbidcat+C2,master)
tabc3.2 <- svytable(~comorbidcat+C3,master)

tabemp0.2 <- svytable(~emp,master)
tabemp1.2 <- svytable(~emp+C1,master)
tabemp2.2 <- svytable(~emp+C2,master)
tabemp3.2 <- svytable(~emp+C3,master)

tabhype0.2 <- svytable(~hypertension,master)
tabhype1.2 <- svytable(~hypertension+C1,master)
tabhype2.2 <- svytable(~hypertension+C2,master)
tabhype3.2 <- svytable(~hypertension+C3,master)

tabchd0.2 <- svytable(~chd,master)
tabchd1.2 <- svytable(~chd+C1,master)
tabchd2.2 <- svytable(~chd+C2,master)
tabchd3.2 <- svytable(~chd+C3,master)

tabang0.2 <- svytable(~angina,master)
tabang1.2 <- svytable(~angina+C1,master)
tabang2.2 <- svytable(~angina+C2,master)
tabang3.2 <- svytable(~angina+C3,master)

tabha0.2 <- svytable(~heartattack,master)
tabha1.2 <- svytable(~heartattack+C1,master)
tabha2.2 <- svytable(~heartattack+C2,master)
tabha3.2 <- svytable(~heartattack+C3,master)

tabhd0.2 <- svytable(~heartdisease,master)
tabhd1.2 <- svytable(~heartdisease+C1,master)
tabhd2.2 <- svytable(~heartdisease+C2,master)
tabhd3.2 <- svytable(~heartdisease+C3,master)

tabstroke0.2 <- svytable(~stroke,master)
tabstroke1.2 <- svytable(~stroke+C1,master)
tabstroke2.2 <- svytable(~stroke+C2,master)
tabstroke3.2 <- svytable(~stroke+C3,master)

tabdiab0.2 <- svytable(~diab,master)
tabdiab1.2 <- svytable(~diab+C1,master)
tabdiab2.2 <- svytable(~diab+C2,master)
tabdiab3.2 <- svytable(~diab+C3,master)

tabbron0.2 <- svytable(~bron,master)
tabbron1.2 <- svytable(~bron+C1,master)
tabbron2.2 <- svytable(~bron+C2,master)
tabbron3.2 <- svytable(~bron+C3,master)

tabkid0.2 <- svytable(~kidney,master)
tabkid1.2 <- svytable(~kidney+C1,master)
tabkid2.2 <- svytable(~kidney+C2,master)
tabkid3.2 <- svytable(~kidney+C3,master)

tabliv0.2 <- svytable(~liver,master)
tabliv1.2 <- svytable(~liver+C1,master)
tabliv2.2 <- svytable(~liver+C2,master)
tabliv3.2 <- svytable(~liver+C3,master)

tabpc0.2 <- svytable(~prior.cancer,master)
tabpc1.2 <- svytable(~prior.cancer+C1,master)
tabpc2.2 <- svytable(~prior.cancer+C2,master)
tabpc3.2 <- svytable(~prior.cancer+C3,master)

tabse0.2 <- svytable(~speceq,master)
tabse1.2 <- svytable(~speceq+C1,master)
tabse2.2 <- svytable(~speceq+C2,master)
tabse3.2 <- svytable(~speceq+C3,master)

agec2.2 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.2 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

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
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.3 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_lcdrat <- subset(master,uspstf.eligible&lcdrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcdrat.eligible&lyg.eligible)
lcdrat_lyg <- subset(master,!uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcdrat.eligible&!lyg.eligible)
lcdrat_only <- subset(master,!uspstf.eligible&lcdrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcdrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)

#Output
#Differences in Selection
tab0.3 <- svytable(~uspstf.eligible|lcdrat.eligible|lyg.eligible,master)
tab1.3 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.3 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.3 <- svytable(~lcdrat.eligible+lyg.eligible,master)
tab4.3 <- svytable(~uspstf.eligible&lcdrat.eligible&lyg.eligible,master)
tab5.3 <- svytable(~uspstf.eligible&lcdrat.eligible&!lyg.eligible,master)
tab6.3 <- svytable(~uspstf.eligible&!lcdrat.eligible&lyg.eligible,master)
tab7.3 <- svytable(~uspstf.eligible==0&lcdrat.eligible&lyg.eligible,master)
tab8.3 <- svytable(~uspstf.eligible&!lcdrat.eligible&!lyg.eligible,master)
tab9.3 <- svytable(~uspstf.eligible==0&lcdrat.eligible&!lyg.eligible,master)
tab10.3 <- svytable(~uspstf.eligible==0&!lcdrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.3 <- svytotal(~LCRAT,master)
lc1.3 <- svytotal(~LCRAT,uspstf)
lc2.3 <- svytotal(~LCRAT,lcdrat)
lc3.3 <- svytotal(~LCRAT,lyg)
lc4.3 <- svytotal(~LCRAT,all)
lc5.3 <- svytotal(~LCRAT,uspstf_lcdrat)
lc6.3 <- svytotal(~LCRAT,uspstf_lyg)
lc7.3 <- svytotal(~LCRAT,lcdrat_lyg)
lc8.3 <- svytotal(~LCRAT,uspstf_only)
lc9.3 <- svytotal(~LCRAT,lcdrat_only)
lc10.3 <- svytotal(~LCRAT,lyg_only)
lc11.3 <- svytotal(~LCRAT,lcdrat_notlyg)
lc12.3 <- svytotal(~LCRAT,lyg_notlcdrat)

#Number of Lung Cancer Deaths
lcd0.3 <- svytotal(~lcd5,master)
lcd1.3 <- svytotal(~lcd5,uspstf)
lcd2.3 <- svytotal(~lcd5,lcdrat)
lcd3.3 <- svytotal(~lcd5,lyg)
lcd4.3 <- svytotal(~lcd5,all)
lcd5.3 <- svytotal(~lcd5,uspstf_lcdrat)
lcd6.3 <- svytotal(~lcd5,uspstf_lyg)
lcd7.3 <- svytotal(~lcd5,lcdrat_lyg)
lcd8.3 <- svytotal(~lcd5,uspstf_only)
lcd9.3 <- svytotal(~lcd5,lcdrat_only)
lcd10.3 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.3 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.3 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.3 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.3 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.3 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.3 <- svytotal(~lcd5,uspstf_lcdrat)*(1-0.796)
lcds6.3 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.3 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds8.3 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.3 <- svytotal(~lcd5,lcdrat_only)*(1-0.796)
lcds10.3 <- svytotal(~lcd5,lyg_only)*(1-0.796)
lcds11.3 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds12.3 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)

#Number of Life-years Gained
lyg0.3 <- svytotal(~lyg,master)
lyg1.3 <- svytotal(~lyg,uspstf)
lyg2.3 <- svytotal(~lyg,lcdrat)
lyg3.3 <- svytotal(~lyg,lyg)
lyg4.3 <- svytotal(~lyg,all)
lyg5.3 <- svytotal(~lyg,uspstf_lcdrat)
lyg6.3 <- svytotal(~lyg,uspstf_lyg)
lyg7.3 <- svytotal(~lyg,lcdrat_lyg)
lyg8.3 <- svytotal(~lyg,uspstf_only)
lyg9.3 <- svytotal(~lyg,lcdrat_only)
lyg10.3 <- svytotal(~lyg,lyg_only)
lyg11.3 <- svytotal(~lyg,lcdrat_notlyg)
lyg12.3 <- svytotal(~lyg,lyg_notlcdrat)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCDRAT.q.3 <- svyquantile(~lcd5,any,quantiles)
nhis$LCDRATcat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
                        ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[2] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[3],2,
                               ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[3] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[4],3,
                                      ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[4] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[5],4,5))))

lyg.q.3 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcdrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCDRAT0.3 <- svymean(~lcd5,master) 
mLCDRAT1.3 <- svymean(~lcd5,c1)
mLCDRAT2.3 <- svymean(~lcd5,c2)
mLCDRAT3.3 <- svymean(~lcd5,c3)
mlyg0.3 <- svymean(~lyg,master) 
mlyg1.3 <- svymean(~lyg,c1)
mlyg2.3 <- svymean(~lyg,c2)
mlyg3.3 <- svymean(~lyg,c3)

tabLCDRAT0.3 <- svytable(~LCDRATcat,master)
tabLCDRAT1.3 <- svytable(~LCDRATcat+C1,master)
tabLCDRAT2.3 <- svytable(~LCDRATcat+C2,master)
tabLCDRAT3.3 <- svytable(~LCDRATcat+C3,master)

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
meda.LCDRAT.3 <- svyquantile(~age,lcdrat,.5)
meda.lyg.3 <- svyquantile(~age,lyg,.5)
meana1.3 <- svymean(~age,c1)
meana2.3 <- svymean(~age,c2)
meana3.3 <- svymean(~age,c3)

taba0.3 <- svytable(~age.cat,master)
taba1.3 <- svytable(~age.cat+C1,master)
taba2.3 <- svytable(~age.cat+C2,master)
taba3.3 <- svytable(~age.cat+C3,master)

tabs0.3 <- svytable(~current,master)
tabs1.3 <- svytable(~current+C1,master)
tabs2.3 <- svytable(~current+C2,master)
tabs3.3 <- svytable(~current+C3,master)

meanc0.3 <- svymean(~comorbidities,master)
meanc1.3 <- svymean(~comorbidities,c1)
meanc2.3 <- svymean(~comorbidities,c2)
meanc3.3 <- svymean(~comorbidities,c3)

tabc0.3 <- svytable(~comorbidcat,master)
tabc1.3 <- svytable(~comorbidcat+C1,master)
tabc2.3 <- svytable(~comorbidcat+C2,master)
tabc3.3 <- svytable(~comorbidcat+C3,master)

tabemp0.3 <- svytable(~emp,master)
tabemp1.3 <- svytable(~emp+C1,master)
tabemp2.3 <- svytable(~emp+C2,master)
tabemp3.3 <- svytable(~emp+C3,master)

tabhype0.3 <- svytable(~hypertension,master)
tabhype1.3 <- svytable(~hypertension+C1,master)
tabhype2.3 <- svytable(~hypertension+C2,master)
tabhype3.3 <- svytable(~hypertension+C3,master)

tabchd0.3 <- svytable(~chd,master)
tabchd1.3 <- svytable(~chd+C1,master)
tabchd2.3 <- svytable(~chd+C2,master)
tabchd3.3 <- svytable(~chd+C3,master)

tabang0.3 <- svytable(~angina,master)
tabang1.3 <- svytable(~angina+C1,master)
tabang2.3 <- svytable(~angina+C2,master)
tabang3.3 <- svytable(~angina+C3,master)

tabha0.3 <- svytable(~heartattack,master)
tabha1.3 <- svytable(~heartattack+C1,master)
tabha2.3 <- svytable(~heartattack+C2,master)
tabha3.3 <- svytable(~heartattack+C3,master)

tabhd0.3 <- svytable(~heartdisease,master)
tabhd1.3 <- svytable(~heartdisease+C1,master)
tabhd2.3 <- svytable(~heartdisease+C2,master)
tabhd3.3 <- svytable(~heartdisease+C3,master)

tabstroke0.3 <- svytable(~stroke,master)
tabstroke1.3 <- svytable(~stroke+C1,master)
tabstroke2.3 <- svytable(~stroke+C2,master)
tabstroke3.3 <- svytable(~stroke+C3,master)

tabdiab0.3 <- svytable(~diab,master)
tabdiab1.3 <- svytable(~diab+C1,master)
tabdiab2.3 <- svytable(~diab+C2,master)
tabdiab3.3 <- svytable(~diab+C3,master)

tabbron0.3 <- svytable(~bron,master)
tabbron1.3 <- svytable(~bron+C1,master)
tabbron2.3 <- svytable(~bron+C2,master)
tabbron3.3 <- svytable(~bron+C3,master)

tabkid0.3 <- svytable(~kidney,master)
tabkid1.3 <- svytable(~kidney+C1,master)
tabkid2.3 <- svytable(~kidney+C2,master)
tabkid3.3 <- svytable(~kidney+C3,master)

tabliv0.3 <- svytable(~liver,master)
tabliv1.3 <- svytable(~liver+C1,master)
tabliv2.3 <- svytable(~liver+C2,master)
tabliv3.3 <- svytable(~liver+C3,master)

tabpc0.3 <- svytable(~prior.cancer,master)
tabpc1.3 <- svytable(~prior.cancer+C1,master)
tabpc2.3 <- svytable(~prior.cancer+C2,master)
tabpc3.3 <- svytable(~prior.cancer+C3,master)

tabse0.3 <- svytable(~speceq,master)
tabse1.3 <- svytable(~speceq+C1,master)
tabse2.3 <- svytable(~speceq+C2,master)
tabse3.3 <- svytable(~speceq+C3,master)

agec2.3 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.3 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

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
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.4 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_lcdrat <- subset(master,uspstf.eligible&lcdrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcdrat.eligible&lyg.eligible)
lcdrat_lyg <- subset(master,!uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcdrat.eligible&!lyg.eligible)
lcdrat_only <- subset(master,!uspstf.eligible&lcdrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcdrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)

#Output
#Differences in Selection
tab0.4 <- svytable(~uspstf.eligible|lcdrat.eligible|lyg.eligible,master)
tab1.4 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.4 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.4 <- svytable(~lcdrat.eligible+lyg.eligible,master)
tab4.4 <- svytable(~uspstf.eligible&lcdrat.eligible&lyg.eligible,master)
tab5.4 <- svytable(~uspstf.eligible&lcdrat.eligible&!lyg.eligible,master)
tab6.4 <- svytable(~uspstf.eligible&!lcdrat.eligible&lyg.eligible,master)
tab7.4 <- svytable(~uspstf.eligible==0&lcdrat.eligible&lyg.eligible,master)
tab8.4 <- svytable(~uspstf.eligible&!lcdrat.eligible&!lyg.eligible,master)
tab9.4 <- svytable(~uspstf.eligible==0&lcdrat.eligible&!lyg.eligible,master)
tab10.4 <- svytable(~uspstf.eligible==0&!lcdrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.4 <- svytotal(~LCRAT,master)
lc1.4 <- svytotal(~LCRAT,uspstf)
lc2.4 <- svytotal(~LCRAT,lcdrat)
lc3.4 <- svytotal(~LCRAT,lyg)
lc4.4 <- svytotal(~LCRAT,all)
lc5.4 <- svytotal(~LCRAT,uspstf_lcdrat)
lc6.4 <- svytotal(~LCRAT,uspstf_lyg)
lc7.4 <- svytotal(~LCRAT,lcdrat_lyg)
lc8.4 <- svytotal(~LCRAT,uspstf_only)
lc9.4 <- svytotal(~LCRAT,lcdrat_only)
lc10.4 <- svytotal(~LCRAT,lyg_only)
lc11.4 <- svytotal(~LCRAT,lcdrat_notlyg)
lc12.4 <- svytotal(~LCRAT,lyg_notlcdrat)

#Number of Lung Cancer Deaths
lcd0.4 <- svytotal(~lcd5,master)
lcd1.4 <- svytotal(~lcd5,uspstf)
lcd2.4 <- svytotal(~lcd5,lcdrat)
lcd3.4 <- svytotal(~lcd5,lyg)
lcd4.4 <- svytotal(~lcd5,all)
lcd5.4 <- svytotal(~lcd5,uspstf_lcdrat)
lcd6.4 <- svytotal(~lcd5,uspstf_lyg)
lcd7.4 <- svytotal(~lcd5,lcdrat_lyg)
lcd8.4 <- svytotal(~lcd5,uspstf_only)
lcd9.4 <- svytotal(~lcd5,lcdrat_only)
lcd10.4 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.4 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.4 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.4 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.4 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.4 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.4 <- svytotal(~lcd5,uspstf_lcdrat)*(1-0.796)
lcds6.4 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.4 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds8.4 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.4 <- svytotal(~lcd5,lcdrat_only)*(1-0.796)
lcds10.4 <- svytotal(~lcd5,lyg_only)*(1-0.796)
lcds11.4 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds12.4 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)

#Number of Life-years Gained
lyg0.4 <- svytotal(~lyg,master)
lyg1.4 <- svytotal(~lyg,uspstf)
lyg2.4 <- svytotal(~lyg,lcdrat)
lyg3.4 <- svytotal(~lyg,lyg)
lyg4.4 <- svytotal(~lyg,all)
lyg5.4 <- svytotal(~lyg,uspstf_lcdrat)
lyg6.4 <- svytotal(~lyg,uspstf_lyg)
lyg7.4 <- svytotal(~lyg,lcdrat_lyg)
lyg8.4 <- svytotal(~lyg,uspstf_only)
lyg9.4 <- svytotal(~lyg,lcdrat_only)
lyg10.4 <- svytotal(~lyg,lyg_only)
lyg11.4 <- svytotal(~lyg,lcdrat_notlyg)
lyg12.4 <- svytotal(~lyg,lyg_notlcdrat)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCDRAT.q.4 <- svyquantile(~lcd5,any,quantiles)
nhis$LCDRATcat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
                        ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[2] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[3],2,
                               ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[3] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[4],3,
                                      ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[4] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[5],4,5))))

lyg.q.4 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcdrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCDRAT0.4 <- svymean(~lcd5,master) 
mLCDRAT1.4 <- svymean(~lcd5,c1)
mLCDRAT2.4 <- svymean(~lcd5,c2)
mLCDRAT3.4 <- svymean(~lcd5,c3)
mlyg0.4 <- svymean(~lyg,master) 
mlyg1.4 <- svymean(~lyg,c1)
mlyg2.4 <- svymean(~lyg,c2)
mlyg3.4 <- svymean(~lyg,c3)

tabLCDRAT0.4 <- svytable(~LCDRATcat,master)
tabLCDRAT1.4 <- svytable(~LCDRATcat+C1,master)
tabLCDRAT2.4 <- svytable(~LCDRATcat+C2,master)
tabLCDRAT3.4 <- svytable(~LCDRATcat+C3,master)

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
meda.LCDRAT.4 <- svyquantile(~age,lcdrat,.5)
meda.lyg.4 <- svyquantile(~age,lyg,.5)
meana1.4 <- svymean(~age,c1)
meana2.4 <- svymean(~age,c2)
meana3.4 <- svymean(~age,c3)

taba0.4 <- svytable(~age.cat,master)
taba1.4 <- svytable(~age.cat+C1,master)
taba2.4 <- svytable(~age.cat+C2,master)
taba3.4 <- svytable(~age.cat+C3,master)

tabs0.4 <- svytable(~current,master)
tabs1.4 <- svytable(~current+C1,master)
tabs2.4 <- svytable(~current+C2,master)
tabs3.4 <- svytable(~current+C3,master)

meanc0.4 <- svymean(~comorbidities,master)
meanc1.4 <- svymean(~comorbidities,c1)
meanc2.4 <- svymean(~comorbidities,c2)
meanc3.4 <- svymean(~comorbidities,c3)

tabc0.4 <- svytable(~comorbidcat,master)
tabc1.4 <- svytable(~comorbidcat+C1,master)
tabc2.4 <- svytable(~comorbidcat+C2,master)
tabc3.4 <- svytable(~comorbidcat+C3,master)

tabemp0.4 <- svytable(~emp,master)
tabemp1.4 <- svytable(~emp+C1,master)
tabemp2.4 <- svytable(~emp+C2,master)
tabemp3.4 <- svytable(~emp+C3,master)

tabhype0.4 <- svytable(~hypertension,master)
tabhype1.4 <- svytable(~hypertension+C1,master)
tabhype2.4 <- svytable(~hypertension+C2,master)
tabhype3.4 <- svytable(~hypertension+C3,master)

tabchd0.4 <- svytable(~chd,master)
tabchd1.4 <- svytable(~chd+C1,master)
tabchd2.4 <- svytable(~chd+C2,master)
tabchd3.4 <- svytable(~chd+C3,master)

tabang0.4 <- svytable(~angina,master)
tabang1.4 <- svytable(~angina+C1,master)
tabang2.4 <- svytable(~angina+C2,master)
tabang3.4 <- svytable(~angina+C3,master)

tabha0.4 <- svytable(~heartattack,master)
tabha1.4 <- svytable(~heartattack+C1,master)
tabha2.4 <- svytable(~heartattack+C2,master)
tabha3.4 <- svytable(~heartattack+C3,master)

tabhd0.4 <- svytable(~heartdisease,master)
tabhd1.4 <- svytable(~heartdisease+C1,master)
tabhd2.4 <- svytable(~heartdisease+C2,master)
tabhd3.4 <- svytable(~heartdisease+C3,master)

tabstroke0.4 <- svytable(~stroke,master)
tabstroke1.4 <- svytable(~stroke+C1,master)
tabstroke2.4 <- svytable(~stroke+C2,master)
tabstroke3.4 <- svytable(~stroke+C3,master)

tabdiab0.4 <- svytable(~diab,master)
tabdiab1.4 <- svytable(~diab+C1,master)
tabdiab2.4 <- svytable(~diab+C2,master)
tabdiab3.4 <- svytable(~diab+C3,master)

tabbron0.4 <- svytable(~bron,master)
tabbron1.4 <- svytable(~bron+C1,master)
tabbron2.4 <- svytable(~bron+C2,master)
tabbron3.4 <- svytable(~bron+C3,master)

tabkid0.4 <- svytable(~kidney,master)
tabkid1.4 <- svytable(~kidney+C1,master)
tabkid2.4 <- svytable(~kidney+C2,master)
tabkid3.4 <- svytable(~kidney+C3,master)

tabliv0.4 <- svytable(~liver,master)
tabliv1.4 <- svytable(~liver+C1,master)
tabliv2.4 <- svytable(~liver+C2,master)
tabliv3.4 <- svytable(~liver+C3,master)

tabpc0.4 <- svytable(~prior.cancer,master)
tabpc1.4 <- svytable(~prior.cancer+C1,master)
tabpc2.4 <- svytable(~prior.cancer+C2,master)
tabpc3.4 <- svytable(~prior.cancer+C3,master)

tabse0.4 <- svytable(~speceq,master)
tabse1.4 <- svytable(~speceq+C1,master)
tabse2.4 <- svytable(~speceq+C2,master)
tabse3.4 <- svytable(~speceq+C3,master)

agec2.4 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.4 <- svytable(~age+C3,master)[,2]/svytable(~age,master)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

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
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000

nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.5 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$lyg.eligible==1,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
all <- subset(master,uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_lcdrat <- subset(master,uspstf.eligible&lcdrat.eligible&!lyg.eligible)
uspstf_lyg <- subset(master,uspstf.eligible&!lcdrat.eligible&lyg.eligible)
lcdrat_lyg <- subset(master,!uspstf.eligible&lcdrat.eligible&lyg.eligible)
uspstf_only <- subset(master,uspstf.eligible&!lcdrat.eligible&!lyg.eligible)
lcdrat_only <- subset(master,!uspstf.eligible&lcdrat.eligible&!lyg.eligible)
lyg_only <- subset(master,!uspstf.eligible&!lcdrat.eligible&lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)

#Output
#Differences in Selection
tab0.5 <- svytable(~uspstf.eligible|lcdrat.eligible|lyg.eligible,master)
tab1.5 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.5 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.5 <- svytable(~lcdrat.eligible+lyg.eligible,master)
tab4.5 <- svytable(~uspstf.eligible&lcdrat.eligible&lyg.eligible,master)
tab5.5 <- svytable(~uspstf.eligible&lcdrat.eligible&!lyg.eligible,master)
tab6.5 <- svytable(~uspstf.eligible&!lcdrat.eligible&lyg.eligible,master)
tab7.5 <- svytable(~uspstf.eligible==0&lcdrat.eligible&lyg.eligible,master)
tab8.5 <- svytable(~uspstf.eligible&!lcdrat.eligible&!lyg.eligible,master)
tab9.5 <- svytable(~uspstf.eligible==0&lcdrat.eligible&!lyg.eligible,master)
tab10.5 <- svytable(~uspstf.eligible==0&!lcdrat.eligible&lyg.eligible,master)

#Number of Lung Cancers
lc0.5 <- svytotal(~LCRAT,master)
lc1.5 <- svytotal(~LCRAT,uspstf)
lc2.5 <- svytotal(~LCRAT,lcdrat)
lc3.5 <- svytotal(~LCRAT,lyg)
lc4.5 <- svytotal(~LCRAT,all)
lc5.5 <- svytotal(~LCRAT,uspstf_lcdrat)
lc6.5 <- svytotal(~LCRAT,uspstf_lyg)
lc7.5 <- svytotal(~LCRAT,lcdrat_lyg)
lc8.5 <- svytotal(~LCRAT,uspstf_only)
lc9.5 <- svytotal(~LCRAT,lcdrat_only)
lc10.5 <- svytotal(~LCRAT,lyg_only)
lc11.5 <- svytotal(~LCRAT,lcdrat_notlyg)
lc12.5 <- svytotal(~LCRAT,lyg_notlcdrat)

#Number of Lung Cancer Deaths
lcd0.5 <- svytotal(~lcd5,master)
lcd1.5 <- svytotal(~lcd5,uspstf)
lcd2.5 <- svytotal(~lcd5,lcdrat)
lcd3.5 <- svytotal(~lcd5,lyg)
lcd4.5 <- svytotal(~lcd5,all)
lcd5.5 <- svytotal(~lcd5,uspstf_lcdrat)
lcd6.5 <- svytotal(~lcd5,uspstf_lyg)
lcd7.5 <- svytotal(~lcd5,lcdrat_lyg)
lcd8.5 <- svytotal(~lcd5,uspstf_only)
lcd9.5 <- svytotal(~lcd5,lcdrat_only)
lcd10.5 <- svytotal(~lcd5,lyg_only)

#Number of Lung Cancer Deaths Saved
lcds0.5 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.5 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.5 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.5 <- svytotal(~lcd5,lyg)*(1-0.796)
lcds4.5 <- svytotal(~lcd5,all)*(1-0.796)
lcds5.5 <- svytotal(~lcd5,uspstf_lcdrat)*(1-0.796)
lcds6.5 <- svytotal(~lcd5,uspstf_lyg)*(1-0.796)
lcds7.5 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds8.5 <- svytotal(~lcd5,uspstf_only)*(1-0.796)
lcds9.5 <- svytotal(~lcd5,lcdrat_only)*(1-0.796)
lcds10.5 <- svytotal(~lcd5,lyg_only)*(1-0.796)
lcds11.5 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds12.5 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)

#Number of Life-years Gained
lyg0.5 <- svytotal(~lyg,master)
lyg1.5 <- svytotal(~lyg,uspstf)
lyg2.5 <- svytotal(~lyg,lcdrat)
lyg3.5 <- svytotal(~lyg,lyg)
lyg4.5 <- svytotal(~lyg,all)
lyg5.5 <- svytotal(~lyg,uspstf_lcdrat)
lyg6.5 <- svytotal(~lyg,uspstf_lyg)
lyg7.5 <- svytotal(~lyg,lcdrat_lyg)
lyg8.5 <- svytotal(~lyg,uspstf_only)
lyg9.5 <- svytotal(~lyg,lcdrat_only)
lyg10.5 <- svytotal(~lyg,lyg_only)
lyg11.5 <- svytotal(~lyg,lcdrat_notlyg)
lyg12.5 <- svytotal(~lyg,lyg_notlcdrat)

quantiles <- c(0,0.2,0.4,0.6,0.8,1)
LCDRAT.q.5 <- svyquantile(~lcd5,any,quantiles)
nhis$LCDRATcat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
                        ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[2] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[3],2,
                               ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[3] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[4],3,
                                      ifelse(nhis$lcd5>svyquantile(~lcd5,any,quantiles)[4] & nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[5],4,5))))

lyg.q.5 <- svyquantile(~lyg,any,quantiles)
nhis$lygcat <- ifelse(nhis$lyg<=svyquantile(~lyg,any,quantiles)[2],1,
                      ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[2] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[3],2,
                             ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[3] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[4],3,
                                    ifelse(nhis$lyg>svyquantile(~lyg,any,quantiles)[4] & nhis$lyg<=svyquantile(~lyg,any,quantiles)[5],4,5))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
c1 <- subset(master,uspstf.eligible==1)
c2 <- subset(master,lcdrat.eligible==1)
c3 <- subset(master,lyg.eligible==1)

mLCDRAT0.5 <- svymean(~lcd5,master) 
mLCDRAT1.5 <- svymean(~lcd5,c1)
mLCDRAT2.5 <- svymean(~lcd5,c2)
mLCDRAT3.5 <- svymean(~lcd5,c3)
mlyg0.5 <- svymean(~lyg,master) 
mlyg1.5 <- svymean(~lyg,c1)
mlyg2.5 <- svymean(~lyg,c2)
mlyg3.5 <- svymean(~lyg,c3)

tabLCDRAT0.5 <- svytable(~LCDRATcat,master)
tabLCDRAT1.5 <- svytable(~LCDRATcat+C1,master)
tabLCDRAT2.5 <- svytable(~LCDRATcat+C2,master)
tabLCDRAT3.5 <- svytable(~LCDRATcat+C3,master)

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
meda.LCDRAT.5 <- svyquantile(~age,lcdrat,.5)
meda.lyg.5 <- svyquantile(~age,lyg,.5)
meana1.5 <- svymean(~age,c1)
meana2.5 <- svymean(~age,c2)
meana3.5 <- svymean(~age,c3)

taba0.5 <- svytable(~age.cat,master)
taba1.5 <- svytable(~age.cat+C1,master)
taba2.5 <- svytable(~age.cat+C2,master)
taba3.5 <- svytable(~age.cat+C3,master)

tabs0.5 <- svytable(~current,master)
tabs1.5 <- svytable(~current+C1,master)
tabs2.5 <- svytable(~current+C2,master)
tabs3.5 <- svytable(~current+C3,master)

meanc0.5 <- svymean(~comorbidities,master)
meanc1.5 <- svymean(~comorbidities,c1)
meanc2.5 <- svymean(~comorbidities,c2)
meanc3.5 <- svymean(~comorbidities,c3)

tabc0.5 <- svytable(~comorbidcat,master)
tabc1.5 <- svytable(~comorbidcat+C1,master)
tabc2.5 <- svytable(~comorbidcat+C2,master)
tabc3.5 <- svytable(~comorbidcat+C3,master)

tabemp0.5 <- svytable(~emp,master)
tabemp1.5 <- svytable(~emp+C1,master)
tabemp2.5 <- svytable(~emp+C2,master)
tabemp3.5 <- svytable(~emp+C3,master)

tabhype0.5 <- svytable(~hypertension,master)
tabhype1.5 <- svytable(~hypertension+C1,master)
tabhype2.5 <- svytable(~hypertension+C2,master)
tabhype3.5 <- svytable(~hypertension+C3,master)

tabchd0.5 <- svytable(~chd,master)
tabchd1.5 <- svytable(~chd+C1,master)
tabchd2.5 <- svytable(~chd+C2,master)
tabchd3.5 <- svytable(~chd+C3,master)

tabang0.5 <- svytable(~angina,master)
tabang1.5 <- svytable(~angina+C1,master)
tabang2.5 <- svytable(~angina+C2,master)
tabang3.5 <- svytable(~angina+C3,master)

tabha0.5 <- svytable(~heartattack,master)
tabha1.5 <- svytable(~heartattack+C1,master)
tabha2.5 <- svytable(~heartattack+C2,master)
tabha3.5 <- svytable(~heartattack+C3,master)

tabhd0.5 <- svytable(~heartdisease,master)
tabhd1.5 <- svytable(~heartdisease+C1,master)
tabhd2.5 <- svytable(~heartdisease+C2,master)
tabhd3.5 <- svytable(~heartdisease+C3,master)

tabstroke0.5 <- svytable(~stroke,master)
tabstroke1.5 <- svytable(~stroke+C1,master)
tabstroke2.5 <- svytable(~stroke+C2,master)
tabstroke3.5 <- svytable(~stroke+C3,master)

tabdiab0.5 <- svytable(~diab,master)
tabdiab1.5 <- svytable(~diab+C1,master)
tabdiab2.5 <- svytable(~diab+C2,master)
tabdiab3.5 <- svytable(~diab+C3,master)

tabbron0.5 <- svytable(~bron,master)
tabbron1.5 <- svytable(~bron+C1,master)
tabbron2.5 <- svytable(~bron+C2,master)
tabbron3.5 <- svytable(~bron+C3,master)

tabkid0.5 <- svytable(~kidney,master)
tabkid1.5 <- svytable(~kidney+C1,master)
tabkid2.5 <- svytable(~kidney+C2,master)
tabkid3.5 <- svytable(~kidney+C3,master)

tabliv0.5 <- svytable(~liver,master)
tabliv1.5 <- svytable(~liver+C1,master)
tabliv2.5 <- svytable(~liver+C2,master)
tabliv3.5 <- svytable(~liver+C3,master)

tabpc0.5 <- svytable(~prior.cancer,master)
tabpc1.5 <- svytable(~prior.cancer+C1,master)
tabpc2.5 <- svytable(~prior.cancer+C2,master)
tabpc3.5 <- svytable(~prior.cancer+C3,master)

tabse0.5 <- svytable(~speceq,master)
tabse1.5 <- svytable(~speceq+C1,master)
tabse2.5 <- svytable(~speceq+C2,master)
tabse3.5 <- svytable(~speceq+C3,master)

agec2.5 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.5 <- svytable(~age+C3,master)[,2]/svytable(~age,master)


#Average values
(LCDRAT.cutoff.1+LCDRAT.cutoff.2+LCDRAT.cutoff.3+LCDRAT.cutoff.4+LCDRAT.cutoff.5)/5
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
(lc0.1+lc0.2+lc0.3+lc0.4+lc0.5)/(5*697074)  #all lung cancers
(lc1.1+lc1.2+lc1.3+lc1.4+lc1.5)/(5*697074)  #uspstf
(lc2.1+lc2.2+lc2.3+lc2.4+lc2.5)/(5*697074)  #lcdrat
(lc3.1+lc3.2+lc3.3+lc3.4+lc3.5)/(5*697074)  #lyg
(lc4.1+lc4.2+lc4.3+lc4.4+lc4.5)/(5*697074)  #all of uspstf, lcdrat, or lyg
(lc5.1+lc5.2+lc5.3+lc5.4+lc5.5)/(5*697074)  #uspstf and lcdrat
(lc6.1+lc6.2+lc6.3+lc6.4+lc6.5)/(5*697074)  #uspstf and lyg
(lc7.1+lc7.2+lc7.3+lc7.4+lc7.5)/(5*697074)  #lcdrat and lyg
(lc8.1+lc8.2+lc8.3+lc8.4+lc8.5)/(5*697074)  #uspstf only
(lc9.1+lc9.2+lc9.3+lc9.4+lc9.5)/(5*697074)  #lcdrat only
(lc10.1+lc10.2+lc10.3+lc10.4+lc10.5)/(5*697074)  #lyg only
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
(lc11.1+lc11.2+lc11.3+lc11.4+lc11.5)/5
(lc12.1+lc12.2+lc12.3+lc12.4+lc12.5)/5
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
(lcds11.1+lcds11.2+lcds11.3+lcds11.4+lcds11.5)/5
(lcds12.1+lcds12.2+lcds12.3+lcds12.4+lcds12.5)/5
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
(lyg11.1+lyg11.2+lyg11.3+lyg11.4+lyg11.5)/5
(lyg12.1+lyg12.2+lyg12.3+lyg12.4+lyg12.5)/5
agec2 <- (agec2.1+agec2.2+agec2.3+agec2.4+agec2.5)/5
agec3 <- (agec3.1+agec3.2+agec3.3+agec3.4+agec3.5)/5


#overall
100*(LCDRAT.q.1+LCDRAT.q.2+LCDRAT.q.3+LCDRAT.q.4+LCDRAT.q.5)/5
365.25*(lyg.q.1+lyg.q.2+lyg.q.3+lyg.q.4+lyg.q.5)/5
(mLCDRAT0.1+mLCDRAT0.2+mLCDRAT0.3+mLCDRAT0.4+mLCDRAT0.5)/5
365.25*(mlyg0.1+mlyg0.2+mlyg0.3+mlyg0.4+mlyg0.5)/5
(tabLCDRAT0.1+tabLCDRAT0.2+tabLCDRAT0.3+tabLCDRAT0.4+tabLCDRAT0.5)/5
(tablyg0.1+tablyg0.2+tablyg0.3+tablyg0.4+tablyg0.5)/5
(tabg0.1+tabg0.2+tabg0.3+tabg0.4+tabg0.5)/5
(tabr0.1+tabr0.2+tabr0.3+tabr0.4+tabr0.5)/5
(meda0.1+meda0.2+meda0.3+meda0.4+meda0.5)/5 #median age
(meda.uspstf.1+meda.uspstf.2+meda.uspstf.3+meda.uspstf.4+meda.uspstf.5)/5
(meda.LCDRAT.1+meda.LCDRAT.2+meda.LCDRAT.3+meda.LCDRAT.4+meda.LCDRAT.5)/5
(meda.lyg.1+meda.lyg.2+meda.lyg.3+meda.lyg.4+meda.lyg.5)/5
(meana1.1+meana1.2+meana1.3+meana1.4+meana1.5)/5
(meana2.1+meana2.2+meana2.3+meana2.4+meana2.5)/5
(meana3.1+meana3.2+meana3.3+meana3.4+meana3.5)/5
(taba0.1+taba0.2+taba0.3+taba0.4+taba0.5)/5 #age categories
(tabc0.1+tabc0.2+tabc0.3+tabc0.4+tabc0.5)/5
(meanc1.1+meanc1.2+meanc1.3+meanc1.4+meanc1.5)/5
(meanc2.1+meanc2.2+meanc2.3+meanc2.4+meanc2.5)/5
(meanc3.1+meanc3.2+meanc3.3+meanc3.4+meanc3.5)/5

#Revised to match new 6 columns
tab4 <- rbind(100*c((mLCDRAT1.1+mLCDRAT1.2+mLCDRAT1.3+mLCDRAT1.4+mLCDRAT1.5)/5,
(mLCDRAT2.1+mLCDRAT2.2+mLCDRAT2.3+mLCDRAT2.4+mLCDRAT2.5)/5,
(mLCDRAT3.1+mLCDRAT3.2+mLCDRAT3.3+mLCDRAT3.4+mLCDRAT3.5)/5),
c(365.25*(mlyg1.1+mlyg1.2+mlyg1.3+mlyg1.4+mlyg1.5)/5,
365.25*(mlyg2.1+mlyg2.2+mlyg2.3+mlyg2.4+mlyg2.5)/5,
365.25*(mlyg3.1+mlyg3.2+mlyg3.3+mlyg3.4+mlyg3.5)/5),
cbind(((tabLCDRAT1.1+tabLCDRAT1.2+tabLCDRAT1.3+tabLCDRAT1.4+tabLCDRAT1.5)/5)[,2], #LCRAT5 overall
((tabLCDRAT2.1+tabLCDRAT2.2+tabLCDRAT2.3+tabLCDRAT2.4+tabLCDRAT2.5)/5)[,2],
((tabLCDRAT3.1+tabLCDRAT3.2+tabLCDRAT3.3+tabLCDRAT3.4+tabLCDRAT3.5)/5)[,2])/1000,
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
cbind(((tabs1.1+tabs1.2+tabs1.3+tabs1.4+tabs1.5)/5)[,2], #smoking status
      ((tabs2.1+tabs2.2+tabs2.3+tabs2.4+tabs2.5)/5)[,2],
      ((tabs3.1+tabs3.2+tabs3.3+tabs3.4+tabs3.5)/5)[,2])/1000,
cbind(((tabc1.1+tabc1.2+tabc1.3+tabc1.4+tabc1.5)/5)[,2], #comorbidities
((tabc2.1+tabc2.2+tabc2.3+tabc2.4+tabc2.5)/5)[,2],
((tabc3.1+tabc3.2+tabc3.3+tabc3.4+tabc3.5)/5)[,2])/1000,
cbind(((tabemp1.1+tabemp1.2+tabemp1.3+tabemp1.4+tabemp1.5)/5)[,2], #emp
      ((tabemp2.1+tabemp2.2+tabemp2.3+tabemp2.4+tabemp2.5)/5)[,2],
      ((tabemp3.1+tabemp3.2+tabemp3.3+tabemp3.4+tabemp3.5)/5)[,2])/1000,
cbind(((tabhype1.1+tabhype1.2+tabhype1.3+tabhype1.4+tabhype1.5)/5)[,2], #hypertension
      ((tabhype2.1+tabhype2.2+tabhype2.3+tabhype2.4+tabhype2.5)/5)[,2],
      ((tabhype3.1+tabhype3.2+tabhype3.3+tabhype3.4+tabhype3.5)/5)[,2])/1000,
cbind(((tabchd1.1+tabchd1.2+tabchd1.3+tabchd1.4+tabchd1.5)/5)[,2], #chronic heart disease
      ((tabchd2.1+tabchd2.2+tabchd2.3+tabchd2.4+tabchd2.5)/5)[,2],
      ((tabchd3.1+tabchd3.2+tabchd3.3+tabchd3.4+tabchd3.5)/5)[,2])/1000,
cbind(((tabang1.1+tabang1.2+tabang1.3+tabang1.4+tabang1.5)/5)[,2], #angina
      ((tabang2.1+tabang2.2+tabang2.3+tabang2.4+tabang2.5)/5)[,2],
      ((tabang3.1+tabang3.2+tabang3.3+tabang3.4+tabang3.5)/5)[,2])/1000,
cbind(((tabha1.1+tabha1.2+tabha1.3+tabha1.4+tabha1.5)/5)[,2], #heart attack
      ((tabha2.1+tabha2.2+tabha2.3+tabha2.4+tabha2.5)/5)[,2],
      ((tabha3.1+tabha3.2+tabha3.3+tabha3.4+tabha3.5)/5)[,2])/1000,
cbind(((tabhd1.1+tabhd1.2+tabhd1.3+tabhd1.4+tabhd1.5)/5)[,2], #heart disease
      ((tabhd2.1+tabhd2.2+tabhd2.3+tabhd2.4+tabhd2.5)/5)[,2],
      ((tabhd3.1+tabhd3.2+tabhd3.3+tabhd3.4+tabhd3.5)/5)[,2])/1000,
cbind(((tabstroke1.1+tabstroke1.2+tabstroke1.3+tabstroke1.4+tabstroke1.5)/5)[,2], #stroke
      ((tabstroke2.1+tabstroke2.2+tabstroke2.3+tabstroke2.4+tabstroke2.5)/5)[,2],
      ((tabstroke3.1+tabstroke3.2+tabstroke3.3+tabstroke3.4+tabstroke3.5)/5)[,2])/1000,
cbind(((tabdiab1.1+tabdiab1.2+tabdiab1.3+tabdiab1.4+tabdiab1.5)/5)[,2], #diabetes
      ((tabdiab2.1+tabdiab2.2+tabdiab2.3+tabdiab2.4+tabdiab2.5)/5)[,2],
      ((tabdiab3.1+tabdiab3.2+tabdiab3.3+tabdiab3.4+tabdiab3.5)/5)[,2])/1000,
cbind(((tabbron1.1+tabbron1.2+tabbron1.3+tabbron1.4+tabbron1.5)/5)[,2], #bronchitis
      ((tabbron2.1+tabbron2.2+tabbron2.3+tabbron2.4+tabbron2.5)/5)[,2],
      ((tabbron3.1+tabbron3.2+tabbron3.3+tabbron3.4+tabbron3.5)/5)[,2])/1000,
cbind(((tabkid1.1+tabkid1.2+tabkid1.3+tabkid1.4+tabkid1.5)/5)[,2], #kidney disease
      ((tabkid2.1+tabkid2.2+tabkid2.3+tabkid2.4+tabkid2.5)/5)[,2],
      ((tabkid3.1+tabkid3.2+tabkid3.3+tabkid3.4+tabkid3.5)/5)[,2])/1000,
cbind(((tabliv1.1+tabliv1.2+tabliv1.3+tabliv1.4+tabliv1.5)/5)[,2], #liver disease
      ((tabliv2.1+tabliv2.2+tabliv2.3+tabliv2.4+tabliv2.5)/5)[,2],
      ((tabliv3.1+tabliv3.2+tabliv3.3+tabliv3.4+tabliv3.5)/5)[,2])/1000,
cbind(((tabpc1.1+tabpc1.2+tabpc1.3+tabpc1.4+tabpc1.5)/5)[,2], #prior cancer
      ((tabpc2.1+tabpc2.2+tabpc2.3+tabpc2.4+tabpc2.5)/5)[,2],
      ((tabpc3.1+tabpc3.2+tabpc3.3+tabpc3.4+tabpc3.5)/5)[,2])/1000,
cbind(((tabse1.1+tabse1.2+tabse1.3+tabse1.4+tabse1.5)/5)[,2], #spec equipment
      ((tabse2.1+tabse2.2+tabse2.3+tabse2.4+tabse2.5)/5)[,2],
      ((tabse3.1+tabse3.2+tabse3.3+tabse3.4+tabse3.5)/5)[,2])/1000)

colnames(tab4) <- c("USPSTF", "LCDRAT", "LYG")
rownames(tab4) <- c("Mean LCDRAT","Mean LDG", 
                    "Q1 LCDRAT", "Q2 LCDRAT", "Q3 LCDRAT", "Q4 LCDRAT", "Q5 LCDRAT",
                    "Q1 LDG", "Q2 LDG", "Q3 LDG", "Q4 LDG", "Q5 LDG",
                    "Males","Females",
                    "Whites","Blacks","Hispanics","Other",
                    "40-49","50-59","60-69","70-79","80-84",
                    "former","current",
                    "No comorbidities","1 comorbidities","2 comordibities","3+ comordibities",
                    "no emp","emp",
                    "no hypertension","hypertension",
                    "no chronic heart disease","chronic heart disease",
                    "no angina","angina",
                    "no heart attack","heart attack",
                    "no heart disease","heart disease",
                    "no stroke","stroke",
                    "no diabetes","diabetes",
                    "no bronchitis","bronchitis",
                    "no kidney disease","kidney disease",
                    "no liver disease","liver disease",
                    "no prior cancer","prior cancer",
                    "no spec equipment","spec equipment")

write.table(round(tab4,2),"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/screen_selections_margins.v7.csv",sep=",")