#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

inf_eo1 <- function(eo1,eo2,eo3,eo4,eo5){
  ave_eo <- (eo1+eo2+eo3+eo4+eo5)/5  #average estimate
  vars.bar <- (eo1[2,]+eo2[2,]+eo3[2,]+eo4[2,]+eo5[2,])/5  #within imputation variance
  v.impute <- ((eo1[1,]-ave_eo[1,])^2+(eo2[1,]-ave_eo[1,])^2+(eo3[1,]-ave_eo[1,])^2+(eo4[1,]-ave_eo[1,])^2+(eo5[1,]-ave_eo[1,])^2)/4  #between-imputation variance
  var.ave_eo <- (vars.bar+(1+1/5)*v.impute)  #variance for average estimate
  #confidence intervals for average estimate
  lower <- ave_eo[1,]-1.96*sqrt(var.ave_eo)
  upper <- ave_eo[1,]+1.96*sqrt(var.ave_eo)
  ave_eo <- rbind(ave_eo[1,],var.ave_eo,lower,upper)
  rownames(ave_eo) <- c("Ave","Var","LCL","UCL")
  return(ave_eo)
}


inf_eo2 <- function(eo1,eo2,eo3,eo4,eo5){
  ave_eo <- (eo1+eo2+eo3+eo4+eo5)/5  #average proportion
  vars.bar <- (eo1[2]+eo2[2]+eo3[2]+eo4[2]+eo5[2])/5  #within imputation variance
  v.impute <- ((eo1[1]-ave_eo[1])^2+(eo2[1]-ave_eo[1])^2+(eo3[1]-ave_eo[1])^2+(eo4[1]-ave_eo[1])^2+(eo5[1]-ave_eo[1])^2)/4  #between-imputation variance
  var.ave_eo <- (vars.bar+(1+1/5)*v.impute)  #variance for average proportion
  logv <- var.ave_eo/ave_eo[1]^2  #variance of log of average proportion
  #confidence intervals for average proportion
  lower <- exp(log(ave_eo[1])-1.96*sqrt(logv))
  upper <- exp(log(ave_eo[1])+1.96*sqrt(logv))
  ave_eo <- rbind(ave_eo[1],var.ave_eo,lower,upper)
  rownames(ave_eo) <- c("Ave","Var","LCL","UCL")
  return(ave_eo)
}

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

nhis$C1 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C3 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C4 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==0,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

nhis$pkyr.cat <- ifelse(nhis$packyears<10,1,
                        ifelse(nhis$packyears<20,2,
                               ifelse(nhis$packyears<30,3,
                                      ifelse(nhis$packyears<40,4,
                                             ifelse(nhis$packyears<50,5,6)))))
nhis$qtyears.cat <- ifelse(nhis$qtyears==0,0,
                           ifelse(nhis$qtyears<=5,1,
                                  ifelse(nhis$qtyears<=10,2,
                                         ifelse(nhis$qtyears<=15,3,
                                                ifelse(nhis$qtyears<=20,4,5)))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

lcdrat_lyg <- subset(master,lcdrat.eligible&lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)
notlcdrat_notlyg <- subset(master,!lcdrat.eligible&!lyg.eligible)

#Output
#Differences in Selection
tab.1 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancers
lc0.1 <- svytotal(~LCRAT,master)
lc1.1 <- svytotal(~LCRAT,lcdrat_lyg)
lc2.1 <- svytotal(~LCRAT,lcdrat_notlyg)
lc3.1 <- svytotal(~LCRAT,lyg_notlcdrat)
lc4.1 <- svytotal(~LCRAT,notlcdrat_notlyg)

#Number of Lung Cancer Deaths
lcd0.1 <- svytotal(~lcd5,master)
lcd1.1 <- svytotal(~lcd5,lcdrat_lyg)
lcd2.1 <- svytotal(~lcd5,lcdrat_notlyg)
lcd3.1 <- svytotal(~lcd5,lyg_notlcdrat)
lcd4.1 <- svytotal(~lcd5,notlcdrat_notlyg)

#Number of Lung Cancer Deaths Saved
lcds0.1 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.1 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds2.1 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds3.1 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)
lcds4.1 <- svytotal(~lcd5,notlcdrat_notlyg)*(1-0.796)

#Number of Life-years Gained
lyg0.1 <- svytotal(~lyg,master)
lyg1.1 <- svytotal(~lyg,lcdrat_lyg)
lyg2.1 <- svytotal(~lyg,lcdrat_notlyg)
lyg3.1 <- svytotal(~lyg,lyg_notlcdrat)
lyg4.1 <- svytotal(~lyg,notlcdrat_notlyg)

mLCDRAT0.1 <- svymean(~lcd5,master) 
mLCDRAT1.1 <- svymean(~lcd5,lcdrat_lyg)
mLCDRAT2.1 <- svymean(~lcd5,lcdrat_notlyg)
mLCDRAT3.1 <- svymean(~lcd5,lyg_notlcdrat)
mLCDRAT4.1 <- svymean(~lcd5,notlcdrat_notlyg)
mlyg0.1 <- svymean(~lyg,master) 
mlyg1.1 <- svymean(~lyg,lcdrat_lyg)
mlyg2.1 <- svymean(~lyg,lcdrat_notlyg)
mlyg3.1 <- svymean(~lyg,lyg_notlcdrat)
mlyg4.1 <- svymean(~lyg,notlcdrat_notlyg)

tabg0.1 <- svytable(~female,master)
tabg1.1 <- svytable(~female+C1,master)
tabg2.1 <- svytable(~female+C2,master)
tabg3.1 <- svytable(~female+C3,master)
tabg4.1 <- svytable(~female+C4,master)

tabr0.1 <- svytable(~race,master)
tabr1.1 <- svytable(~race+C1,master)
tabr2.1 <- svytable(~race+C2,master)
tabr3.1 <- svytable(~race+C3,master)
tabr4.1 <- svytable(~race+C4,master)

meana0.1 <- svymean(~age,lcdrat_lyg)
meana1.1 <- svymean(~age,lcdrat_lyg)
meana2.1 <- svymean(~age,lcdrat_notlyg)
meana3.1 <- svymean(~age,lyg_notlcdrat)
meana4.1 <- svymean(~age,notlcdrat_notlyg)

taba0.1 <- svytable(~age.cat,master)
taba1.1 <- svytable(~age.cat+C1,master)
taba2.1 <- svytable(~age.cat+C2,master)
taba3.1 <- svytable(~age.cat+C3,master)
taba4.1 <- svytable(~age.cat+C4,master)

tabs0.1 <- svytable(~current,master)
tabs1.1 <- svytable(~current+C1,master)
tabs2.1 <- svytable(~current+C2,master)
tabs3.1 <- svytable(~current+C3,master)
tabs4.1 <- svytable(~current+C4,master)

tabpk0.1 <- svytable(~pkyr.cat,master)
tabpk1.1 <- svytable(~pkyr.cat+C1,master)
tabpk2.1 <- svytable(~pkyr.cat+C2,master)
tabpk3.1 <- svytable(~pkyr.cat+C3,master)
tabpk4.1 <- svytable(~pkyr.cat+C4,master)

tabqt0.1 <- svytable(~qtyears.cat,master)
tabqt1.1 <- svytable(~qtyears.cat+C1,master)
tabqt2.1 <- svytable(~qtyears.cat+C2,master)
tabqt3.1 <- svytable(~qtyears.cat+C3,master)
tabqt4.1 <- svytable(~qtyears.cat+C4,master)

meanc0.1 <- svymean(~comorbidities,master)
meanc1.1 <- svymean(~comorbidities,lcdrat_lyg)
meanc2.1 <- svymean(~comorbidities,lcdrat_notlyg)
meanc3.1 <- svymean(~comorbidities,lyg_notlcdrat)
meanc4.1 <- svymean(~comorbidities,notlcdrat_notlyg)

tabc0.1 <- svytable(~comorbidcat,master)
tabc1.1 <- svytable(~comorbidcat+C1,master)
tabc2.1 <- svytable(~comorbidcat+C2,master)
tabc3.1 <- svytable(~comorbidcat+C3,master)
tabc4.1 <- svytable(~comorbidcat+C4,master)

tabemp0.1 <- svytable(~emp,master)
tabemp1.1 <- svytable(~emp+C1,master)
tabemp2.1 <- svytable(~emp+C2,master)
tabemp3.1 <- svytable(~emp+C3,master)
tabemp4.1 <- svytable(~emp+C4,master)

tabhype0.1 <- svytable(~hypertension,master)
tabhype1.1 <- svytable(~hypertension+C1,master)
tabhype2.1 <- svytable(~hypertension+C2,master)
tabhype3.1 <- svytable(~hypertension+C3,master)
tabhype4.1 <- svytable(~hypertension+C4,master)

tabchd0.1 <- svytable(~chd,master)
tabchd1.1 <- svytable(~chd+C1,master)
tabchd2.1 <- svytable(~chd+C2,master)
tabchd3.1 <- svytable(~chd+C3,master)
tabchd4.1 <- svytable(~chd+C4,master)

tabang0.1 <- svytable(~angina,master)
tabang1.1 <- svytable(~angina+C1,master)
tabang2.1 <- svytable(~angina+C2,master)
tabang3.1 <- svytable(~angina+C3,master)
tabang4.1 <- svytable(~angina+C4,master)

tabha0.1 <- svytable(~heartattack,master)
tabha1.1 <- svytable(~heartattack+C1,master)
tabha2.1 <- svytable(~heartattack+C2,master)
tabha3.1 <- svytable(~heartattack+C3,master)
tabha4.1 <- svytable(~heartattack+C4,master)

tabhd0.1 <- svytable(~heartdisease,master)
tabhd1.1 <- svytable(~heartdisease+C1,master)
tabhd2.1 <- svytable(~heartdisease+C2,master)
tabhd3.1 <- svytable(~heartdisease+C3,master)
tabhd4.1 <- svytable(~heartdisease+C4,master)

tabstroke0.1 <- svytable(~stroke,master)
tabstroke1.1 <- svytable(~stroke+C1,master)
tabstroke2.1 <- svytable(~stroke+C2,master)
tabstroke3.1 <- svytable(~stroke+C3,master)
tabstroke4.1 <- svytable(~stroke+C4,master)

tabdiab0.1 <- svytable(~diab,master)
tabdiab1.1 <- svytable(~diab+C1,master)
tabdiab2.1 <- svytable(~diab+C2,master)
tabdiab3.1 <- svytable(~diab+C3,master)
tabdiab4.1 <- svytable(~diab+C4,master)

tabbron0.1 <- svytable(~bron,master)
tabbron1.1 <- svytable(~bron+C1,master)
tabbron2.1 <- svytable(~bron+C2,master)
tabbron3.1 <- svytable(~bron+C3,master)
tabbron4.1 <- svytable(~bron+C4,master)

tabkid0.1 <- svytable(~kidney,master)
tabkid1.1 <- svytable(~kidney+C1,master)
tabkid2.1 <- svytable(~kidney+C2,master)
tabkid3.1 <- svytable(~kidney+C3,master)
tabkid4.1 <- svytable(~kidney+C4,master)

tabliv0.1 <- svytable(~liver,master)
tabliv1.1 <- svytable(~liver+C1,master)
tabliv2.1 <- svytable(~liver+C2,master)
tabliv3.1 <- svytable(~liver+C3,master)
tabliv4.1 <- svytable(~liver+C4,master)

tabpc0.1 <- svytable(~prior.cancer,master)
tabpc1.1 <- svytable(~prior.cancer+C1,master)
tabpc2.1 <- svytable(~prior.cancer+C2,master)
tabpc3.1 <- svytable(~prior.cancer+C3,master)
tabpc4.1 <- svytable(~prior.cancer+C4,master)

tabse0.1 <- svytable(~speceq,master)
tabse1.1 <- svytable(~speceq+C1,master)
tabse2.1 <- svytable(~speceq+C2,master)
tabse3.1 <- svytable(~speceq+C3,master)
tabse4.1 <- svytable(~speceq+C4,master)

agec1.1 <- svytable(~age+C1,master)[,2]/svytable(~age,master)
agec2.1 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.1 <- svytable(~age+C3,master)[,2]/svytable(~age,master)
agec4.1 <- svytable(~age+C4,master)[,2]/svytable(~age,master)

NNS1.1 <- svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)
NNS2.1 <- svyratio(~lyg.eligible,~lcds5,lyg_notlcdrat)
NNS3.1 <- svyratio(~lcdrat.eligible,~lyg,lcdrat_notlyg)
NNS4.1 <- svyratio(~lyg.eligible,~lyg,lyg_notlcdrat)
LGLCDS1.1 <- svyratio(~lyg,~lcds5,lcdrat_notlyg)
LGLCDS2.1 <- svyratio(~lyg,~lcds5,lyg_notlcdrat)

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


nhis$C1 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C3 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C4 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==0,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

nhis$pkyr.cat <- ifelse(nhis$packyears<10,1,
                        ifelse(nhis$packyears<20,2,
                               ifelse(nhis$packyears<30,3,
                                      ifelse(nhis$packyears<40,4,
                                             ifelse(nhis$packyears<50,5,6)))))
nhis$qtyears.cat <- ifelse(nhis$qtyears==0,0,
                           ifelse(nhis$qtyears<=5,1,
                                  ifelse(nhis$qtyears<=10,2,
                                         ifelse(nhis$qtyears<=15,3,
                                                ifelse(nhis$qtyears<=20,4,5)))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

lcdrat_lyg <- subset(master,lcdrat.eligible&lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)
notlcdrat_notlyg <- subset(master,!lcdrat.eligible&!lyg.eligible)

#Output
#Differences in Selection
tab.2 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancers
lc0.2 <- svytotal(~LCRAT,master)
lc1.2 <- svytotal(~LCRAT,lcdrat_lyg)
lc2.2 <- svytotal(~LCRAT,lcdrat_notlyg)
lc3.2 <- svytotal(~LCRAT,lyg_notlcdrat)
lc4.2 <- svytotal(~LCRAT,notlcdrat_notlyg)

#Number of Lung Cancer Deaths
lcd0.2 <- svytotal(~lcd5,master)
lcd1.2 <- svytotal(~lcd5,lcdrat_lyg)
lcd2.2 <- svytotal(~lcd5,lcdrat_notlyg)
lcd3.2 <- svytotal(~lcd5,lyg_notlcdrat)
lcd4.2 <- svytotal(~lcd5,notlcdrat_notlyg)

#Number of Lung Cancer Deaths Saved
lcds0.2 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.2 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds2.2 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds3.2 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)
lcds4.2 <- svytotal(~lcd5,notlcdrat_notlyg)*(1-0.796)

#Number of Life-years Gained
lyg0.2 <- svytotal(~lyg,master)
lyg1.2 <- svytotal(~lyg,lcdrat_lyg)
lyg2.2 <- svytotal(~lyg,lcdrat_notlyg)
lyg3.2 <- svytotal(~lyg,lyg_notlcdrat)
lyg4.2 <- svytotal(~lyg,notlcdrat_notlyg)

mLCDRAT0.2 <- svymean(~lcd5,master) 
mLCDRAT1.2 <- svymean(~lcd5,lcdrat_lyg)
mLCDRAT2.2 <- svymean(~lcd5,lcdrat_notlyg)
mLCDRAT3.2 <- svymean(~lcd5,lyg_notlcdrat)
mLCDRAT4.2 <- svymean(~lcd5,notlcdrat_notlyg)
mlyg0.2 <- svymean(~lyg,master) 
mlyg1.2 <- svymean(~lyg,lcdrat_lyg)
mlyg2.2 <- svymean(~lyg,lcdrat_notlyg)
mlyg3.2 <- svymean(~lyg,lyg_notlcdrat)
mlyg4.2 <- svymean(~lyg,notlcdrat_notlyg)

svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)


tabg0.2 <- svytable(~female,master)
tabg1.2 <- svytable(~female+C1,master)
tabg2.2 <- svytable(~female+C2,master)
tabg3.2 <- svytable(~female+C3,master)
tabg4.2 <- svytable(~female+C4,master)

tabr0.2 <- svytable(~race,master)
tabr1.2 <- svytable(~race+C1,master)
tabr2.2 <- svytable(~race+C2,master)
tabr3.2 <- svytable(~race+C3,master)
tabr4.2 <- svytable(~race+C4,master)

meana0.2 <- svymean(~age,lcdrat_lyg)
meana1.2 <- svymean(~age,lcdrat_lyg)
meana2.2 <- svymean(~age,lcdrat_notlyg)
meana3.2 <- svymean(~age,lyg_notlcdrat)
meana4.2 <- svymean(~age,notlcdrat_notlyg)

taba0.2 <- svytable(~age.cat,master)
taba1.2 <- svytable(~age.cat+C1,master)
taba2.2 <- svytable(~age.cat+C2,master)
taba3.2 <- svytable(~age.cat+C3,master)
taba4.2 <- svytable(~age.cat+C4,master)

tabs0.2 <- svytable(~current,master)
tabs1.2 <- svytable(~current+C1,master)
tabs2.2 <- svytable(~current+C2,master)
tabs3.2 <- svytable(~current+C3,master)
tabs4.2 <- svytable(~current+C4,master)

tabpk0.2 <- svytable(~pkyr.cat,master)
tabpk1.2 <- svytable(~pkyr.cat+C1,master)
tabpk2.2 <- svytable(~pkyr.cat+C2,master)
tabpk3.2 <- svytable(~pkyr.cat+C3,master)
tabpk4.2 <- svytable(~pkyr.cat+C4,master)

tabqt0.2 <- svytable(~qtyears.cat,master)
tabqt1.2 <- svytable(~qtyears.cat+C1,master)
tabqt2.2 <- svytable(~qtyears.cat+C2,master)
tabqt3.2 <- svytable(~qtyears.cat+C3,master)
tabqt4.2 <- svytable(~qtyears.cat+C4,master)

meanc0.2 <- svymean(~comorbidities,master)
meanc1.2 <- svymean(~comorbidities,lcdrat_lyg)
meanc2.2 <- svymean(~comorbidities,lcdrat_notlyg)
meanc3.2 <- svymean(~comorbidities,lyg_notlcdrat)
meanc4.2 <- svymean(~comorbidities,notlcdrat_notlyg)

tabc0.2 <- svytable(~comorbidcat,master)
tabc1.2 <- svytable(~comorbidcat+C1,master)
tabc2.2 <- svytable(~comorbidcat+C2,master)
tabc3.2 <- svytable(~comorbidcat+C3,master)
tabc4.2 <- svytable(~comorbidcat+C4,master)

tabemp0.2 <- svytable(~emp,master)
tabemp1.2 <- svytable(~emp+C1,master)
tabemp2.2 <- svytable(~emp+C2,master)
tabemp3.2 <- svytable(~emp+C3,master)
tabemp4.2 <- svytable(~emp+C4,master)

tabhype0.2 <- svytable(~hypertension,master)
tabhype1.2 <- svytable(~hypertension+C1,master)
tabhype2.2 <- svytable(~hypertension+C2,master)
tabhype3.2 <- svytable(~hypertension+C3,master)
tabhype4.2 <- svytable(~hypertension+C4,master)

tabchd0.2 <- svytable(~chd,master)
tabchd1.2 <- svytable(~chd+C1,master)
tabchd2.2 <- svytable(~chd+C2,master)
tabchd3.2 <- svytable(~chd+C3,master)
tabchd4.2 <- svytable(~chd+C4,master)

tabang0.2 <- svytable(~angina,master)
tabang1.2 <- svytable(~angina+C1,master)
tabang2.2 <- svytable(~angina+C2,master)
tabang3.2 <- svytable(~angina+C3,master)
tabang4.2 <- svytable(~angina+C4,master)

tabha0.2 <- svytable(~heartattack,master)
tabha1.2 <- svytable(~heartattack+C1,master)
tabha2.2 <- svytable(~heartattack+C2,master)
tabha3.2 <- svytable(~heartattack+C3,master)
tabha4.2 <- svytable(~heartattack+C4,master)

tabhd0.2 <- svytable(~heartdisease,master)
tabhd1.2 <- svytable(~heartdisease+C1,master)
tabhd2.2 <- svytable(~heartdisease+C2,master)
tabhd3.2 <- svytable(~heartdisease+C3,master)
tabhd4.2 <- svytable(~heartdisease+C4,master)

tabstroke0.2 <- svytable(~stroke,master)
tabstroke1.2 <- svytable(~stroke+C1,master)
tabstroke2.2 <- svytable(~stroke+C2,master)
tabstroke3.2 <- svytable(~stroke+C3,master)
tabstroke4.2 <- svytable(~stroke+C4,master)

tabdiab0.2 <- svytable(~diab,master)
tabdiab1.2 <- svytable(~diab+C1,master)
tabdiab2.2 <- svytable(~diab+C2,master)
tabdiab3.2 <- svytable(~diab+C3,master)
tabdiab4.2 <- svytable(~diab+C4,master)

tabbron0.2 <- svytable(~bron,master)
tabbron1.2 <- svytable(~bron+C1,master)
tabbron2.2 <- svytable(~bron+C2,master)
tabbron3.2 <- svytable(~bron+C3,master)
tabbron4.2 <- svytable(~bron+C4,master)

tabkid0.2 <- svytable(~kidney,master)
tabkid1.2 <- svytable(~kidney+C1,master)
tabkid2.2 <- svytable(~kidney+C2,master)
tabkid3.2 <- svytable(~kidney+C3,master)
tabkid4.2 <- svytable(~kidney+C4,master)

tabliv0.2 <- svytable(~liver,master)
tabliv1.2 <- svytable(~liver+C1,master)
tabliv2.2 <- svytable(~liver+C2,master)
tabliv3.2 <- svytable(~liver+C3,master)
tabliv4.2 <- svytable(~liver+C4,master)

tabpc0.2 <- svytable(~prior.cancer,master)
tabpc1.2 <- svytable(~prior.cancer+C1,master)
tabpc2.2 <- svytable(~prior.cancer+C2,master)
tabpc3.2 <- svytable(~prior.cancer+C3,master)
tabpc4.2 <- svytable(~prior.cancer+C4,master)

tabse0.2 <- svytable(~speceq,master)
tabse1.2 <- svytable(~speceq+C1,master)
tabse2.2 <- svytable(~speceq+C2,master)
tabse3.2 <- svytable(~speceq+C3,master)
tabse4.2 <- svytable(~speceq+C4,master)

agec1.2 <- svytable(~age+C1,master)[,2]/svytable(~age,master)
agec2.2 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.2 <- svytable(~age+C3,master)[,2]/svytable(~age,master)
agec4.2 <- svytable(~age+C4,master)[,2]/svytable(~age,master)

NNS1.2 <- svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)
NNS2.2 <- svyratio(~lyg.eligible,~lcds5,lyg_notlcdrat)
NNS3.2 <- svyratio(~lcdrat.eligible,~lyg,lcdrat_notlyg)
NNS4.2 <- svyratio(~lyg.eligible,~lyg,lyg_notlcdrat)
LGLCDS1.2 <- svyratio(~lyg,~lcds5,lcdrat_notlyg)
LGLCDS2.2 <- svyratio(~lyg,~lcds5,lyg_notlcdrat)

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


nhis$C1 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C3 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C4 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==0,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

nhis$pkyr.cat <- ifelse(nhis$packyears<10,1,
                        ifelse(nhis$packyears<20,2,
                               ifelse(nhis$packyears<30,3,
                                      ifelse(nhis$packyears<40,4,
                                             ifelse(nhis$packyears<50,5,6)))))
nhis$qtyears.cat <- ifelse(nhis$qtyears==0,0,
                           ifelse(nhis$qtyears<=5,1,
                                  ifelse(nhis$qtyears<=10,2,
                                         ifelse(nhis$qtyears<=15,3,
                                                ifelse(nhis$qtyears<=20,4,5)))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

lcdrat_lyg <- subset(master,lcdrat.eligible&lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)
notlcdrat_notlyg <- subset(master,!lcdrat.eligible&!lyg.eligible)

#Output
#Differences in Selection
tab.3 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancers
lc0.3 <- svytotal(~LCRAT,master)
lc1.3 <- svytotal(~LCRAT,lcdrat_lyg)
lc2.3 <- svytotal(~LCRAT,lcdrat_notlyg)
lc3.3 <- svytotal(~LCRAT,lyg_notlcdrat)
lc4.3 <- svytotal(~LCRAT,notlcdrat_notlyg)

#Number of Lung Cancer Deaths
lcd0.3 <- svytotal(~lcd5,master)
lcd1.3 <- svytotal(~lcd5,lcdrat_lyg)
lcd2.3 <- svytotal(~lcd5,lcdrat_notlyg)
lcd3.3 <- svytotal(~lcd5,lyg_notlcdrat)
lcd4.3 <- svytotal(~lcd5,notlcdrat_notlyg)

#Number of Lung Cancer Deaths Saved
lcds0.3 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.3 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds2.3 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds3.3 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)
lcds4.3 <- svytotal(~lcd5,notlcdrat_notlyg)*(1-0.796)

#Number of Life-years Gained
lyg0.3 <- svytotal(~lyg,master)
lyg1.3 <- svytotal(~lyg,lcdrat_lyg)
lyg2.3 <- svytotal(~lyg,lcdrat_notlyg)
lyg3.3 <- svytotal(~lyg,lyg_notlcdrat)
lyg4.3 <- svytotal(~lyg,notlcdrat_notlyg)

mLCDRAT0.3 <- svymean(~lcd5,master) 
mLCDRAT1.3 <- svymean(~lcd5,lcdrat_lyg)
mLCDRAT2.3 <- svymean(~lcd5,lcdrat_notlyg)
mLCDRAT3.3 <- svymean(~lcd5,lyg_notlcdrat)
mLCDRAT4.3 <- svymean(~lcd5,notlcdrat_notlyg)
mlyg0.3 <- svymean(~lyg,master) 
mlyg1.3 <- svymean(~lyg,lcdrat_lyg)
mlyg2.3 <- svymean(~lyg,lcdrat_notlyg)
mlyg3.3 <- svymean(~lyg,lyg_notlcdrat)
mlyg4.3 <- svymean(~lyg,notlcdrat_notlyg)

tabg0.3 <- svytable(~female,master)
tabg1.3 <- svytable(~female+C1,master)
tabg2.3 <- svytable(~female+C2,master)
tabg3.3 <- svytable(~female+C3,master)
tabg4.3 <- svytable(~female+C4,master)

tabr0.3 <- svytable(~race,master)
tabr1.3 <- svytable(~race+C1,master)
tabr2.3 <- svytable(~race+C2,master)
tabr3.3 <- svytable(~race+C3,master)
tabr4.3 <- svytable(~race+C4,master)

meana0.3 <- svymean(~age,lcdrat_lyg)
meana1.3 <- svymean(~age,lcdrat_lyg)
meana2.3 <- svymean(~age,lcdrat_notlyg)
meana3.3 <- svymean(~age,lyg_notlcdrat)
meana4.3 <- svymean(~age,notlcdrat_notlyg)

taba0.3 <- svytable(~age.cat,master)
taba1.3 <- svytable(~age.cat+C1,master)
taba2.3 <- svytable(~age.cat+C2,master)
taba3.3 <- svytable(~age.cat+C3,master)
taba4.3 <- svytable(~age.cat+C4,master)

tabs0.3 <- svytable(~current,master)
tabs1.3 <- svytable(~current+C1,master)
tabs2.3 <- svytable(~current+C2,master)
tabs3.3 <- svytable(~current+C3,master)
tabs4.3 <- svytable(~current+C4,master)

tabpk0.3 <- svytable(~pkyr.cat,master)
tabpk1.3 <- svytable(~pkyr.cat+C1,master)
tabpk2.3 <- svytable(~pkyr.cat+C2,master)
tabpk3.3 <- svytable(~pkyr.cat+C3,master)
tabpk4.3 <- svytable(~pkyr.cat+C4,master)

tabqt0.3 <- svytable(~qtyears.cat,master)
tabqt1.3 <- svytable(~qtyears.cat+C1,master)
tabqt2.3 <- svytable(~qtyears.cat+C2,master)
tabqt3.3 <- svytable(~qtyears.cat+C3,master)
tabqt4.3 <- svytable(~qtyears.cat+C4,master)

meanc0.3 <- svymean(~comorbidities,master)
meanc1.3 <- svymean(~comorbidities,lcdrat_lyg)
meanc2.3 <- svymean(~comorbidities,lcdrat_notlyg)
meanc3.3 <- svymean(~comorbidities,lyg_notlcdrat)
meanc4.3 <- svymean(~comorbidities,notlcdrat_notlyg)

tabc0.3 <- svytable(~comorbidcat,master)
tabc1.3 <- svytable(~comorbidcat+C1,master)
tabc2.3 <- svytable(~comorbidcat+C2,master)
tabc3.3 <- svytable(~comorbidcat+C3,master)
tabc4.3 <- svytable(~comorbidcat+C4,master)

tabemp0.3 <- svytable(~emp,master)
tabemp1.3 <- svytable(~emp+C1,master)
tabemp2.3 <- svytable(~emp+C2,master)
tabemp3.3 <- svytable(~emp+C3,master)
tabemp4.3 <- svytable(~emp+C4,master)

tabhype0.3 <- svytable(~hypertension,master)
tabhype1.3 <- svytable(~hypertension+C1,master)
tabhype2.3 <- svytable(~hypertension+C2,master)
tabhype3.3 <- svytable(~hypertension+C3,master)
tabhype4.3 <- svytable(~hypertension+C4,master)

tabchd0.3 <- svytable(~chd,master)
tabchd1.3 <- svytable(~chd+C1,master)
tabchd2.3 <- svytable(~chd+C2,master)
tabchd3.3 <- svytable(~chd+C3,master)
tabchd4.3 <- svytable(~chd+C4,master)

tabang0.3 <- svytable(~angina,master)
tabang1.3 <- svytable(~angina+C1,master)
tabang2.3 <- svytable(~angina+C2,master)
tabang3.3 <- svytable(~angina+C3,master)
tabang4.3 <- svytable(~angina+C4,master)

tabha0.3 <- svytable(~heartattack,master)
tabha1.3 <- svytable(~heartattack+C1,master)
tabha2.3 <- svytable(~heartattack+C2,master)
tabha3.3 <- svytable(~heartattack+C3,master)
tabha4.3 <- svytable(~heartattack+C4,master)

tabhd0.3 <- svytable(~heartdisease,master)
tabhd1.3 <- svytable(~heartdisease+C1,master)
tabhd2.3 <- svytable(~heartdisease+C2,master)
tabhd3.3 <- svytable(~heartdisease+C3,master)
tabhd4.3 <- svytable(~heartdisease+C4,master)

tabstroke0.3 <- svytable(~stroke,master)
tabstroke1.3 <- svytable(~stroke+C1,master)
tabstroke2.3 <- svytable(~stroke+C2,master)
tabstroke3.3 <- svytable(~stroke+C3,master)
tabstroke4.3 <- svytable(~stroke+C4,master)

tabdiab0.3 <- svytable(~diab,master)
tabdiab1.3 <- svytable(~diab+C1,master)
tabdiab2.3 <- svytable(~diab+C2,master)
tabdiab3.3 <- svytable(~diab+C3,master)
tabdiab4.3 <- svytable(~diab+C4,master)

tabbron0.3 <- svytable(~bron,master)
tabbron1.3 <- svytable(~bron+C1,master)
tabbron2.3 <- svytable(~bron+C2,master)
tabbron3.3 <- svytable(~bron+C3,master)
tabbron4.3 <- svytable(~bron+C4,master)

tabkid0.3 <- svytable(~kidney,master)
tabkid1.3 <- svytable(~kidney+C1,master)
tabkid2.3 <- svytable(~kidney+C2,master)
tabkid3.3 <- svytable(~kidney+C3,master)
tabkid4.3 <- svytable(~kidney+C4,master)

tabliv0.3 <- svytable(~liver,master)
tabliv1.3 <- svytable(~liver+C1,master)
tabliv2.3 <- svytable(~liver+C2,master)
tabliv3.3 <- svytable(~liver+C3,master)
tabliv4.3 <- svytable(~liver+C4,master)

tabpc0.3 <- svytable(~prior.cancer,master)
tabpc1.3 <- svytable(~prior.cancer+C1,master)
tabpc2.3 <- svytable(~prior.cancer+C2,master)
tabpc3.3 <- svytable(~prior.cancer+C3,master)
tabpc4.3 <- svytable(~prior.cancer+C4,master)

tabse0.3 <- svytable(~speceq,master)
tabse1.3 <- svytable(~speceq+C1,master)
tabse2.3 <- svytable(~speceq+C2,master)
tabse3.3 <- svytable(~speceq+C3,master)
tabse4.3 <- svytable(~speceq+C4,master)

agec1.3 <- svytable(~age+C1,master)[,2]/svytable(~age,master)
agec2.3 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.3 <- svytable(~age+C3,master)[,2]/svytable(~age,master)
agec4.3 <- svytable(~age+C4,master)[,2]/svytable(~age,master)

NNS1.3 <- svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)
NNS2.3 <- svyratio(~lyg.eligible,~lcds5,lyg_notlcdrat)
NNS3.3 <- svyratio(~lcdrat.eligible,~lyg,lcdrat_notlyg)
NNS4.3 <- svyratio(~lyg.eligible,~lyg,lyg_notlcdrat)
LGLCDS1.3 <- svyratio(~lyg,~lcds5,lcdrat_notlyg)
LGLCDS2.3 <- svyratio(~lyg,~lcds5,lyg_notlcdrat)

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


nhis$C1 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C3 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C4 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==0,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

nhis$pkyr.cat <- ifelse(nhis$packyears<10,1,
                        ifelse(nhis$packyears<20,2,
                               ifelse(nhis$packyears<30,3,
                                      ifelse(nhis$packyears<40,4,
                                             ifelse(nhis$packyears<50,5,6)))))
nhis$qtyears.cat <- ifelse(nhis$qtyears==0,0,
                           ifelse(nhis$qtyears<=5,1,
                                  ifelse(nhis$qtyears<=10,2,
                                         ifelse(nhis$qtyears<=15,3,
                                                ifelse(nhis$qtyears<=20,4,5)))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

lcdrat_lyg <- subset(master,lcdrat.eligible&lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)
notlcdrat_notlyg <- subset(master,!lcdrat.eligible&!lyg.eligible)

#Output
#Differences in Selection
tab.4 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancers
lc0.4 <- svytotal(~LCRAT,master)
lc1.4 <- svytotal(~LCRAT,lcdrat_lyg)
lc2.4 <- svytotal(~LCRAT,lcdrat_notlyg)
lc3.4 <- svytotal(~LCRAT,lyg_notlcdrat)
lc4.4 <- svytotal(~LCRAT,notlcdrat_notlyg)

#Number of Lung Cancer Deaths
lcd0.4 <- svytotal(~lcd5,master)
lcd1.4 <- svytotal(~lcd5,lcdrat_lyg)
lcd2.4 <- svytotal(~lcd5,lcdrat_notlyg)
lcd3.4 <- svytotal(~lcd5,lyg_notlcdrat)
lcd4.4 <- svytotal(~lcd5,notlcdrat_notlyg)

#Number of Lung Cancer Deaths Saved
lcds0.4 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.4 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds2.4 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds3.4 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)
lcds4.4 <- svytotal(~lcd5,notlcdrat_notlyg)*(1-0.796)

#Number of Life-years Gained
lyg0.4 <- svytotal(~lyg,master)
lyg1.4 <- svytotal(~lyg,lcdrat_lyg)
lyg2.4 <- svytotal(~lyg,lcdrat_notlyg)
lyg3.4 <- svytotal(~lyg,lyg_notlcdrat)
lyg4.4 <- svytotal(~lyg,notlcdrat_notlyg)

mLCDRAT0.4 <- svymean(~lcd5,master) 
mLCDRAT1.4 <- svymean(~lcd5,lcdrat_lyg)
mLCDRAT2.4 <- svymean(~lcd5,lcdrat_notlyg)
mLCDRAT3.4 <- svymean(~lcd5,lyg_notlcdrat)
mLCDRAT4.4 <- svymean(~lcd5,notlcdrat_notlyg)
mlyg0.4 <- svymean(~lyg,master) 
mlyg1.4 <- svymean(~lyg,lcdrat_lyg)
mlyg2.4 <- svymean(~lyg,lcdrat_notlyg)
mlyg3.4 <- svymean(~lyg,lyg_notlcdrat)
mlyg4.4 <- svymean(~lyg,notlcdrat_notlyg)

tabg0.4 <- svytable(~female,master)
tabg1.4 <- svytable(~female+C1,master)
tabg2.4 <- svytable(~female+C2,master)
tabg3.4 <- svytable(~female+C3,master)
tabg4.4 <- svytable(~female+C4,master)

tabr0.4 <- svytable(~race,master)
tabr1.4 <- svytable(~race+C1,master)
tabr2.4 <- svytable(~race+C2,master)
tabr3.4 <- svytable(~race+C3,master)
tabr4.4 <- svytable(~race+C4,master)

meana0.4 <- svymean(~age,lcdrat_lyg)
meana1.4 <- svymean(~age,lcdrat_lyg)
meana2.4 <- svymean(~age,lcdrat_notlyg)
meana3.4 <- svymean(~age,lyg_notlcdrat)
meana4.4 <- svymean(~age,notlcdrat_notlyg)

taba0.4 <- svytable(~age.cat,master)
taba1.4 <- svytable(~age.cat+C1,master)
taba2.4 <- svytable(~age.cat+C2,master)
taba3.4 <- svytable(~age.cat+C3,master)
taba4.4 <- svytable(~age.cat+C4,master)

tabs0.4 <- svytable(~current,master)
tabs1.4 <- svytable(~current+C1,master)
tabs2.4 <- svytable(~current+C2,master)
tabs3.4 <- svytable(~current+C3,master)
tabs4.4 <- svytable(~current+C4,master)

tabpk0.4 <- svytable(~pkyr.cat,master)
tabpk1.4 <- svytable(~pkyr.cat+C1,master)
tabpk2.4 <- svytable(~pkyr.cat+C2,master)
tabpk3.4 <- svytable(~pkyr.cat+C3,master)
tabpk4.4 <- svytable(~pkyr.cat+C4,master)

tabqt0.4 <- svytable(~qtyears.cat,master)
tabqt1.4 <- svytable(~qtyears.cat+C1,master)
tabqt2.4 <- svytable(~qtyears.cat+C2,master)
tabqt3.4 <- svytable(~qtyears.cat+C3,master)
tabqt4.4 <- svytable(~qtyears.cat+C4,master)

meanc0.4 <- svymean(~comorbidities,master)
meanc1.4 <- svymean(~comorbidities,lcdrat_lyg)
meanc2.4 <- svymean(~comorbidities,lcdrat_notlyg)
meanc3.4 <- svymean(~comorbidities,lyg_notlcdrat)
meanc4.4 <- svymean(~comorbidities,notlcdrat_notlyg)

tabc0.4 <- svytable(~comorbidcat,master)
tabc1.4 <- svytable(~comorbidcat+C1,master)
tabc2.4 <- svytable(~comorbidcat+C2,master)
tabc3.4 <- svytable(~comorbidcat+C3,master)
tabc4.4 <- svytable(~comorbidcat+C4,master)

tabemp0.4 <- svytable(~emp,master)
tabemp1.4 <- svytable(~emp+C1,master)
tabemp2.4 <- svytable(~emp+C2,master)
tabemp3.4 <- svytable(~emp+C3,master)
tabemp4.4 <- svytable(~emp+C4,master)

tabhype0.4 <- svytable(~hypertension,master)
tabhype1.4 <- svytable(~hypertension+C1,master)
tabhype2.4 <- svytable(~hypertension+C2,master)
tabhype3.4 <- svytable(~hypertension+C3,master)
tabhype4.4 <- svytable(~hypertension+C4,master)

tabchd0.4 <- svytable(~chd,master)
tabchd1.4 <- svytable(~chd+C1,master)
tabchd2.4 <- svytable(~chd+C2,master)
tabchd3.4 <- svytable(~chd+C3,master)
tabchd4.4 <- svytable(~chd+C4,master)

tabang0.4 <- svytable(~angina,master)
tabang1.4 <- svytable(~angina+C1,master)
tabang2.4 <- svytable(~angina+C2,master)
tabang3.4 <- svytable(~angina+C3,master)
tabang4.4 <- svytable(~angina+C4,master)

tabha0.4 <- svytable(~heartattack,master)
tabha1.4 <- svytable(~heartattack+C1,master)
tabha2.4 <- svytable(~heartattack+C2,master)
tabha3.4 <- svytable(~heartattack+C3,master)
tabha4.4 <- svytable(~heartattack+C4,master)

tabhd0.4 <- svytable(~heartdisease,master)
tabhd1.4 <- svytable(~heartdisease+C1,master)
tabhd2.4 <- svytable(~heartdisease+C2,master)
tabhd3.4 <- svytable(~heartdisease+C3,master)
tabhd4.4 <- svytable(~heartdisease+C4,master)

tabstroke0.4 <- svytable(~stroke,master)
tabstroke1.4 <- svytable(~stroke+C1,master)
tabstroke2.4 <- svytable(~stroke+C2,master)
tabstroke3.4 <- svytable(~stroke+C3,master)
tabstroke4.4 <- svytable(~stroke+C4,master)

tabdiab0.4 <- svytable(~diab,master)
tabdiab1.4 <- svytable(~diab+C1,master)
tabdiab2.4 <- svytable(~diab+C2,master)
tabdiab3.4 <- svytable(~diab+C3,master)
tabdiab4.4 <- svytable(~diab+C4,master)

tabbron0.4 <- svytable(~bron,master)
tabbron1.4 <- svytable(~bron+C1,master)
tabbron2.4 <- svytable(~bron+C2,master)
tabbron3.4 <- svytable(~bron+C3,master)
tabbron4.4 <- svytable(~bron+C4,master)

tabkid0.4 <- svytable(~kidney,master)
tabkid1.4 <- svytable(~kidney+C1,master)
tabkid2.4 <- svytable(~kidney+C2,master)
tabkid3.4 <- svytable(~kidney+C3,master)
tabkid4.4 <- svytable(~kidney+C4,master)

tabliv0.4 <- svytable(~liver,master)
tabliv1.4 <- svytable(~liver+C1,master)
tabliv2.4 <- svytable(~liver+C2,master)
tabliv3.4 <- svytable(~liver+C3,master)
tabliv4.4 <- svytable(~liver+C4,master)

tabpc0.4 <- svytable(~prior.cancer,master)
tabpc1.4 <- svytable(~prior.cancer+C1,master)
tabpc2.4 <- svytable(~prior.cancer+C2,master)
tabpc3.4 <- svytable(~prior.cancer+C3,master)
tabpc4.4 <- svytable(~prior.cancer+C4,master)

tabse0.4 <- svytable(~speceq,master)
tabse1.4 <- svytable(~speceq+C1,master)
tabse2.4 <- svytable(~speceq+C2,master)
tabse3.4 <- svytable(~speceq+C3,master)
tabse4.4 <- svytable(~speceq+C4,master)

agec1.4 <- svytable(~age+C1,master)[,2]/svytable(~age,master)
agec2.4 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.4 <- svytable(~age+C3,master)[,2]/svytable(~age,master)
agec4.4 <- svytable(~age+C4,master)[,2]/svytable(~age,master)

NNS1.4 <- svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)
NNS2.4 <- svyratio(~lyg.eligible,~lcds5,lyg_notlcdrat)
NNS3.4 <- svyratio(~lcdrat.eligible,~lyg,lcdrat_notlyg)
NNS4.4 <- svyratio(~lyg.eligible,~lyg,lyg_notlcdrat)
LGLCDS1.4 <- svyratio(~lyg,~lcds5,lcdrat_notlyg)
LGLCDS2.4 <- svyratio(~lyg,~lcds5,lyg_notlcdrat)

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


nhis$C1 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==1,1,0)
nhis$C2 <- ifelse(nhis$lcdrat.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C3 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C4 <- ifelse(nhis$lcdrat.eligible==0 & nhis$lyg.eligible==0,1,0)

nhis$age.cat <- ifelse(nhis$age>=80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))

nhis$comorbidcat <- ifelse(nhis$comorbidities==0,0,
                           ifelse(nhis$comorbidities==1,1,
                                  ifelse(nhis$comorbidities==2,2,3)))

nhis$pkyr.cat <- ifelse(nhis$packyears<10,1,
                        ifelse(nhis$packyears<20,2,
                               ifelse(nhis$packyears<30,3,
                                      ifelse(nhis$packyears<40,4,
                                             ifelse(nhis$packyears<50,5,6)))))
nhis$qtyears.cat <- ifelse(nhis$qtyears==0,0,
                           ifelse(nhis$qtyears<=5,1,
                                  ifelse(nhis$qtyears<=10,2,
                                         ifelse(nhis$qtyears<=15,3,
                                                ifelse(nhis$qtyears<=20,4,5)))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

lcdrat_lyg <- subset(master,lcdrat.eligible&lyg.eligible)
lcdrat_notlyg <- subset(master,lcdrat.eligible&!lyg.eligible)
lyg_notlcdrat <- subset(master,!lcdrat.eligible&lyg.eligible)
notlcdrat_notlyg <- subset(master,!lcdrat.eligible&!lyg.eligible)

#Output
#Differences in Selection
tab.5 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancers
lc0.5 <- svytotal(~LCRAT,master)
lc1.5 <- svytotal(~LCRAT,lcdrat_lyg)
lc2.5 <- svytotal(~LCRAT,lcdrat_notlyg)
lc3.5 <- svytotal(~LCRAT,lyg_notlcdrat)
lc4.5 <- svytotal(~LCRAT,notlcdrat_notlyg)

#Number of Lung Cancer Deaths
lcd0.5 <- svytotal(~lcd5,master)
lcd1.5 <- svytotal(~lcd5,lcdrat_lyg)
lcd2.5 <- svytotal(~lcd5,lcdrat_notlyg)
lcd3.5 <- svytotal(~lcd5,lyg_notlcdrat)
lcd4.5 <- svytotal(~lcd5,notlcdrat_notlyg)

#Number of Lung Cancer Deaths Saved
lcds0.5 <- svytotal(~lcd5,master)*(1-0.796)
lcds1.5 <- svytotal(~lcd5,lcdrat_lyg)*(1-0.796)
lcds2.5 <- svytotal(~lcd5,lcdrat_notlyg)*(1-0.796)
lcds3.5 <- svytotal(~lcd5,lyg_notlcdrat)*(1-0.796)
lcds4.5 <- svytotal(~lcd5,notlcdrat_notlyg)*(1-0.796)

#Number of Life-years Gained
lyg0.5 <- svytotal(~lyg,master)
lyg1.5 <- svytotal(~lyg,lcdrat_lyg)
lyg2.5 <- svytotal(~lyg,lcdrat_notlyg)
lyg3.5 <- svytotal(~lyg,lyg_notlcdrat)
lyg4.5 <- svytotal(~lyg,notlcdrat_notlyg)

mLCDRAT0.5 <- svymean(~lcd5,master) 
mLCDRAT1.5 <- svymean(~lcd5,lcdrat_lyg)
mLCDRAT2.5 <- svymean(~lcd5,lcdrat_notlyg)
mLCDRAT3.5 <- svymean(~lcd5,lyg_notlcdrat)
mLCDRAT4.5 <- svymean(~lcd5,notlcdrat_notlyg)
mlyg0.5 <- svymean(~lyg,master) 
mlyg1.5 <- svymean(~lyg,lcdrat_lyg)
mlyg2.5 <- svymean(~lyg,lcdrat_notlyg)
mlyg3.5 <- svymean(~lyg,lyg_notlcdrat)
mlyg4.5 <- svymean(~lyg,notlcdrat_notlyg)

tabg0.5 <- svytable(~female,master)
tabg1.5 <- svytable(~female+C1,master)
tabg2.5 <- svytable(~female+C2,master)
tabg3.5 <- svytable(~female+C3,master)
tabg4.5 <- svytable(~female+C4,master)

tabr0.5 <- svytable(~race,master)
tabr1.5 <- svytable(~race+C1,master)
tabr2.5 <- svytable(~race+C2,master)
tabr3.5 <- svytable(~race+C3,master)
tabr4.5 <- svytable(~race+C4,master)

meana0.5 <- svymean(~age,lcdrat_lyg)
meana1.5 <- svymean(~age,lcdrat_lyg)
meana2.5 <- svymean(~age,lcdrat_notlyg)
meana3.5 <- svymean(~age,lyg_notlcdrat)
meana4.5 <- svymean(~age,notlcdrat_notlyg)

taba0.5 <- svytable(~age.cat,master)
taba1.5 <- svytable(~age.cat+C1,master)
taba2.5 <- svytable(~age.cat+C2,master)
taba3.5 <- svytable(~age.cat+C3,master)
taba4.5 <- svytable(~age.cat+C4,master)

tabs0.5 <- svytable(~current,master)
tabs1.5 <- svytable(~current+C1,master)
tabs2.5 <- svytable(~current+C2,master)
tabs3.5 <- svytable(~current+C3,master)
tabs4.5 <- svytable(~current+C4,master)

tabpk0.5 <- svytable(~pkyr.cat,master)
tabpk1.5 <- svytable(~pkyr.cat+C1,master)
tabpk2.5 <- svytable(~pkyr.cat+C2,master)
tabpk3.5 <- svytable(~pkyr.cat+C3,master)
tabpk4.5 <- svytable(~pkyr.cat+C4,master)

tabqt0.5 <- svytable(~qtyears.cat,master)
tabqt1.5 <- svytable(~qtyears.cat+C1,master)
tabqt2.5 <- svytable(~qtyears.cat+C2,master)
tabqt3.5 <- svytable(~qtyears.cat+C3,master)
tabqt4.5 <- svytable(~qtyears.cat+C4,master)

meanc0.5 <- svymean(~comorbidities,master)
meanc1.5 <- svymean(~comorbidities,lcdrat_lyg)
meanc2.5 <- svymean(~comorbidities,lcdrat_notlyg)
meanc3.5 <- svymean(~comorbidities,lyg_notlcdrat)
meanc4.5 <- svymean(~comorbidities,notlcdrat_notlyg)

tabc0.5 <- svytable(~comorbidcat,master)
tabc1.5 <- svytable(~comorbidcat+C1,master)
tabc2.5 <- svytable(~comorbidcat+C2,master)
tabc3.5 <- svytable(~comorbidcat+C3,master)
tabc4.5 <- svytable(~comorbidcat+C4,master)

tabemp0.5 <- svytable(~emp,master)
tabemp1.5 <- svytable(~emp+C1,master)
tabemp2.5 <- svytable(~emp+C2,master)
tabemp3.5 <- svytable(~emp+C3,master)
tabemp4.5 <- svytable(~emp+C4,master)

tabhype0.5 <- svytable(~hypertension,master)
tabhype1.5 <- svytable(~hypertension+C1,master)
tabhype2.5 <- svytable(~hypertension+C2,master)
tabhype3.5 <- svytable(~hypertension+C3,master)
tabhype4.5 <- svytable(~hypertension+C4,master)

tabchd0.5 <- svytable(~chd,master)
tabchd1.5 <- svytable(~chd+C1,master)
tabchd2.5 <- svytable(~chd+C2,master)
tabchd3.5 <- svytable(~chd+C3,master)
tabchd4.5 <- svytable(~chd+C4,master)

tabang0.5 <- svytable(~angina,master)
tabang1.5 <- svytable(~angina+C1,master)
tabang2.5 <- svytable(~angina+C2,master)
tabang3.5 <- svytable(~angina+C3,master)
tabang4.5 <- svytable(~angina+C4,master)

tabha0.5 <- svytable(~heartattack,master)
tabha1.5 <- svytable(~heartattack+C1,master)
tabha2.5 <- svytable(~heartattack+C2,master)
tabha3.5 <- svytable(~heartattack+C3,master)
tabha4.5 <- svytable(~heartattack+C4,master)

tabhd0.5 <- svytable(~heartdisease,master)
tabhd1.5 <- svytable(~heartdisease+C1,master)
tabhd2.5 <- svytable(~heartdisease+C2,master)
tabhd3.5 <- svytable(~heartdisease+C3,master)
tabhd4.5 <- svytable(~heartdisease+C4,master)

tabstroke0.5 <- svytable(~stroke,master)
tabstroke1.5 <- svytable(~stroke+C1,master)
tabstroke2.5 <- svytable(~stroke+C2,master)
tabstroke3.5 <- svytable(~stroke+C3,master)
tabstroke4.5 <- svytable(~stroke+C4,master)

tabdiab0.5 <- svytable(~diab,master)
tabdiab1.5 <- svytable(~diab+C1,master)
tabdiab2.5 <- svytable(~diab+C2,master)
tabdiab3.5 <- svytable(~diab+C3,master)
tabdiab4.5 <- svytable(~diab+C4,master)

tabbron0.5 <- svytable(~bron,master)
tabbron1.5 <- svytable(~bron+C1,master)
tabbron2.5 <- svytable(~bron+C2,master)
tabbron3.5 <- svytable(~bron+C3,master)
tabbron4.5 <- svytable(~bron+C4,master)

tabkid0.5 <- svytable(~kidney,master)
tabkid1.5 <- svytable(~kidney+C1,master)
tabkid2.5 <- svytable(~kidney+C2,master)
tabkid3.5 <- svytable(~kidney+C3,master)
tabkid4.5 <- svytable(~kidney+C4,master)

tabliv0.5 <- svytable(~liver,master)
tabliv1.5 <- svytable(~liver+C1,master)
tabliv2.5 <- svytable(~liver+C2,master)
tabliv3.5 <- svytable(~liver+C3,master)
tabliv4.5 <- svytable(~liver+C4,master)

tabpc0.5 <- svytable(~prior.cancer,master)
tabpc1.5 <- svytable(~prior.cancer+C1,master)
tabpc2.5 <- svytable(~prior.cancer+C2,master)
tabpc3.5 <- svytable(~prior.cancer+C3,master)
tabpc4.5 <- svytable(~prior.cancer+C4,master)

tabse0.5 <- svytable(~speceq,master)
tabse1.5 <- svytable(~speceq+C1,master)
tabse2.5 <- svytable(~speceq+C2,master)
tabse3.5 <- svytable(~speceq+C3,master)
tabse4.5 <- svytable(~speceq+C4,master)

agec1.5 <- svytable(~age+C1,master)[,2]/svytable(~age,master)
agec2.5 <- svytable(~age+C2,master)[,2]/svytable(~age,master)
agec3.5 <- svytable(~age+C3,master)[,2]/svytable(~age,master)
agec4.5 <- svytable(~age+C4,master)[,2]/svytable(~age,master)

NNS1.5 <- svyratio(~lcdrat.eligible,~lcds5,lcdrat_notlyg)
NNS2.5 <- svyratio(~lyg.eligible,~lcds5,lyg_notlcdrat)
NNS3.5 <- svyratio(~lcdrat.eligible,~lyg,lcdrat_notlyg)
NNS4.5 <- svyratio(~lyg.eligible,~lyg,lyg_notlcdrat)
LGLCDS1.5 <- svyratio(~lyg,~lcds5,lcdrat_notlyg)
LGLCDS2.5 <- svyratio(~lyg,~lcds5,lyg_notlcdrat)

#Average values
(LCDRAT.cutoff.1+LCDRAT.cutoff.2+LCDRAT.cutoff.3+LCDRAT.cutoff.4+LCDRAT.cutoff.5)/5
365.25*(lyg.cutoff.1+lyg.cutoff.2+lyg.cutoff.3+lyg.cutoff.4+lyg.cutoff.5)/5
(tab.1+tab.2+tab.3+tab.4+tab.5)/5
(lc0.1+lc0.2+lc0.3+lc0.4+lc0.5)/5
(lc1.1+lc1.2+lc1.3+lc1.4+lc1.5)/5
(lc2.1+lc2.2+lc2.3+lc2.4+lc2.5)/5
(lc3.1+lc3.2+lc3.3+lc3.4+lc3.5)/5
(lc4.1+lc4.2+lc4.3+lc4.4+lc4.5)/5
(lcd0.1+lcd0.2+lcd0.3+lcd0.4+lcd0.5)/5
(lcd1.1+lcd1.2+lcd1.3+lcd1.4+lcd1.5)/5
(lcd2.1+lcd2.2+lcd2.3+lcd2.4+lcd2.5)/5
(lcd3.1+lcd3.2+lcd3.3+lcd3.4+lcd3.5)/5
(lcd4.1+lcd4.2+lcd4.3+lcd4.4+lcd4.5)/5
(lcds0.1+lcds0.2+lcds0.3+lcds0.4+lcds0.5)/5
(lcds1.1+lcds1.2+lcds1.3+lcds1.4+lcds1.5)/5
(lcds2.1+lcds2.2+lcds2.3+lcds2.4+lcds2.5)/5
(lcds3.1+lcds3.2+lcds3.3+lcds3.4+lcds3.5)/5
(lcds4.1+lcds4.2+lcds4.3+lcds4.4+lcds4.5)/5
(lyg0.1+lyg0.2+lyg0.3+lyg0.4+lyg0.5)/5
(lyg1.1+lyg1.2+lyg1.3+lyg1.4+lyg1.5)/5
(lyg2.1+lyg2.2+lyg2.3+lyg2.4+lyg2.5)/5
(lyg3.1+lyg3.2+lyg3.3+lyg3.4+lyg3.5)/5
(lyg4.1+lyg4.2+lyg4.3+lyg4.4+lyg4.5)/5
agec1 <- (agec1.1+agec1.2+agec1.3+agec1.4+agec1.5)/5
agec2 <- (agec2.1+agec2.2+agec2.3+agec2.4+agec2.5)/5
agec3 <- (agec3.1+agec3.2+agec3.3+agec3.4+agec3.5)/5
agec4 <- (agec4.1+agec4.2+agec4.3+agec4.4+agec4.5)/5

#overall
(mLCDRAT0.1+mLCDRAT0.2+mLCDRAT0.3+mLCDRAT0.4+mLCDRAT0.5)/5
365.25*(mlyg0.1+mlyg0.2+mlyg0.3+mlyg0.4+mlyg0.5)/5
(tabg0.1+tabg0.2+tabg0.3+tabg0.4+tabg0.5)/5
(tabr0.1+tabr0.2+tabr0.3+tabr0.4+tabr0.5)/5
(meana0.1+meana0.2+meana0.3+meana0.4+meana0.5)/5
(meana1.1+meana1.2+meana1.3+meana1.4+meana1.5)/5
(meana2.1+meana2.2+meana2.3+meana2.4+meana2.5)/5
(meana3.1+meana3.2+meana3.3+meana3.4+meana3.5)/5
(meana4.1+meana4.2+meana4.3+meana4.4+meana4.5)/5
(taba0.1+taba0.2+taba0.3+taba0.4+taba0.5)/5 #age categories
(tabc0.1+tabc0.2+tabc0.3+tabc0.4+tabc0.5)/5
(meanc0.1+meanc0.2+meanc0.3+meanc0.4+meanc0.5)/5
(meanc1.1+meanc1.2+meanc1.3+meanc1.4+meanc1.5)/5
(meanc2.1+meanc2.2+meanc2.3+meanc2.4+meanc2.5)/5
(meanc3.1+meanc3.2+meanc3.3+meanc3.4+meanc3.5)/5

inf_eo2(c(mLCDRAT2.1[1],SE(mLCDRAT2.1)^2),
        c(mLCDRAT2.2[1],SE(mLCDRAT2.2)^2),
        c(mLCDRAT2.3[1],SE(mLCDRAT2.3)^2),
        c(mLCDRAT2.4[1],SE(mLCDRAT2.4)^2),
        c(mLCDRAT2.5[1],SE(mLCDRAT2.5)^2))

inf_eo2(c(mLCDRAT3.1[1],SE(mLCDRAT3.1)^2),
        c(mLCDRAT3.2[1],SE(mLCDRAT3.2)^2),
        c(mLCDRAT3.3[1],SE(mLCDRAT3.3)^2),
        c(mLCDRAT3.4[1],SE(mLCDRAT3.4)^2),
        c(mLCDRAT3.5[1],SE(mLCDRAT3.5)^2))

inf_eo1(rbind(c(NNS1.1[[1]],NNS2.1[[1]],NNS3.1[[1]],NNS4.1[[1]],LGLCDS1.1[[1]],LGLCDS2.1[[1]]),
              c(SE(NNS1.1),SE(NNS2.1),SE(NNS3.1),SE(NNS4.1),SE(LGLCDS1.1),SE(LGLCDS2.1))^2),
        rbind(c(NNS1.2[[1]],NNS2.2[[1]],NNS3.2[[1]],NNS4.2[[1]],LGLCDS1.2[[1]],LGLCDS2.2[[1]]),
              c(SE(NNS1.2),SE(NNS2.2),SE(NNS3.2),SE(NNS4.2),SE(LGLCDS1.2),SE(LGLCDS2.2))^2),
        rbind(c(NNS1.3[[1]],NNS2.3[[1]],NNS3.3[[1]],NNS4.3[[1]],LGLCDS1.3[[1]],LGLCDS2.3[[1]]),
              c(SE(NNS1.3),SE(NNS2.3),SE(NNS3.3),SE(NNS4.3),SE(LGLCDS1.3),SE(LGLCDS2.3))^2),
        rbind(c(NNS1.4[[1]],NNS2.4[[1]],NNS3.4[[1]],NNS4.4[[1]],LGLCDS1.4[[1]],LGLCDS2.4[[1]]),
              c(SE(NNS1.4),SE(NNS2.4),SE(NNS3.4),SE(NNS4.4),SE(LGLCDS1.4),SE(LGLCDS2.4))^2),
        rbind(c(NNS1.5[[1]],NNS2.5[[1]],NNS3.5[[1]],NNS4.5[[1]],LGLCDS1.5[[1]],LGLCDS2.5[[1]]),
              c(SE(NNS1.5),SE(NNS2.5),SE(NNS3.5),SE(NNS4.5),SE(LGLCDS1.5),SE(LGLCDS2.5))^2))

#Revised to match new 6 columns
tab4 <- rbind(100*c((mLCDRAT1.1+mLCDRAT1.2+mLCDRAT1.3+mLCDRAT1.4+mLCDRAT1.5)/5,
                    (mLCDRAT2.1+mLCDRAT2.2+mLCDRAT2.3+mLCDRAT2.4+mLCDRAT2.5)/5,
                    (mLCDRAT3.1+mLCDRAT3.2+mLCDRAT3.3+mLCDRAT3.4+mLCDRAT3.5)/5,
                    (mLCDRAT4.1+mLCDRAT4.2+mLCDRAT4.3+mLCDRAT4.4+mLCDRAT4.5)/5),
              c(365.25*(mlyg1.1+mlyg1.2+mlyg1.3+mlyg1.4+mlyg1.5)/5,
                365.25*(mlyg2.1+mlyg2.2+mlyg2.3+mlyg2.4+mlyg2.5)/5,
                365.25*(mlyg3.1+mlyg3.2+mlyg3.3+mlyg3.4+mlyg3.5)/5,
                365.25*(mlyg4.1+mlyg4.2+mlyg4.3+mlyg4.4+mlyg4.5)/5),
              cbind(((tabg1.1+tabg1.2+tabg1.3+tabg1.4+tabg1.5)/5)[,2],  #gender
                    ((tabg2.1+tabg2.2+tabg2.3+tabg2.4+tabg2.5)/5)[,2],
                    ((tabg3.1+tabg3.2+tabg3.3+tabg3.4+tabg3.5)/5)[,2],
                    ((tabg4.1+tabg4.2+tabg4.3+tabg4.4+tabg4.5)/5)[,2])/1000,
              cbind(((tabr1.1+tabr1.2+tabr1.3+tabr1.4+tabr1.5)/5)[,2],  #race
                    ((tabr2.1+tabr2.2+tabr2.3+tabr2.4+tabr2.5)/5)[,2],
                    ((tabr3.1+tabr3.2+tabr3.3+tabr3.4+tabr3.5)/5)[,2],
                    ((tabr4.1+tabr4.2+tabr4.3+tabr4.4+tabr4.5)/5)[,2])/1000,
              c((meana1.1+meana1.2+meana1.3+meana1.4+meana1.5)/5,
                (meana2.1+meana2.2+meana2.3+meana2.4+meana2.5)/5,
                (meana3.1+meana3.2+meana3.3+meana3.4+meana3.5)/5,
                (meana4.1+meana4.2+meana4.3+meana4.4+meana4.5)/5),              
              cbind(((taba1.1+taba1.2+taba1.3+taba1.4+taba1.5)/5)[,2],
                    ((taba2.1+taba2.2+taba2.3+taba2.4+taba2.5)/5)[,2],
                    ((taba3.1+taba3.2+taba3.3+taba3.4+taba3.5)/5)[,2],
                    ((taba4.1+taba4.2+taba4.3+taba4.4+taba4.5)/5)[,2])/1000,
              cbind(((tabs1.1+tabs1.2+tabs1.3+tabs1.4+tabs1.5)/5)[,2], #smoking status
                    ((tabs2.1+tabs2.2+tabs2.3+tabs2.4+tabs2.5)/5)[,2],
                    ((tabs3.1+tabs3.2+tabs3.3+tabs3.4+tabs3.5)/5)[,2],
                    ((tabs4.1+tabs4.2+tabs4.3+tabs4.4+tabs4.5)/5)[,2])/1000,
              cbind(((tabpk1.1+tabpk1.2+tabpk1.3+tabpk1.4+tabpk1.5)/5)[,2], #pk-years
                    ((tabpk2.1+tabpk2.2+tabpk2.3+tabpk2.4+tabpk2.5)/5)[,2],
                    ((tabpk3.1+tabpk3.2+tabpk3.3+tabpk3.4+tabpk3.5)/5)[,2],
                    ((tabpk4.1+tabpk4.2+tabpk4.3+tabpk4.4+tabpk4.5)/5)[,2])/1000,   
              cbind(((tabqt1.1+tabqt1.2+tabqt1.3+tabqt1.4+tabqt1.5)/5)[,2], #qt-years
                    ((tabqt2.1+tabqt2.2+tabqt2.3+tabqt2.4+tabqt2.5)/5)[,2],
                    ((tabqt3.1+tabqt3.2+tabqt3.3+tabqt3.4+tabqt3.5)/5)[,2],
                    ((tabqt4.1+tabqt4.2+tabqt4.3+tabqt4.4+tabqt4.5)/5)[,2])/1000,        
              c((meanc1.1+meanc1.2+meanc1.3+meanc1.4+meanc1.5)/5,
                (meanc2.1+meanc2.2+meanc2.3+meanc2.4+meanc2.5)/5,
                (meanc3.1+meanc3.2+meanc3.3+meanc3.4+meanc3.5)/5,
                (meanc4.1+meanc4.2+meanc4.3+meanc4.4+meanc4.5)/5),                  
              cbind(((tabc1.1+tabc1.2+tabc1.3+tabc1.4+tabc1.5)/5)[,2], #comorbidities
                    ((tabc2.1+tabc2.2+tabc2.3+tabc2.4+tabc2.5)/5)[,2],
                    ((tabc3.1+tabc3.2+tabc3.3+tabc3.4+tabc3.5)/5)[,2],
                    ((tabc4.1+tabc4.2+tabc4.3+tabc4.4+tabc4.5)/5)[,2])/1000,
              cbind(((tabemp1.1+tabemp1.2+tabemp1.3+tabemp1.4+tabemp1.5)/5)[,2], #emp
                    ((tabemp2.1+tabemp2.2+tabemp2.3+tabemp2.4+tabemp2.5)/5)[,2],
                    ((tabemp3.1+tabemp3.2+tabemp3.3+tabemp3.4+tabemp3.5)/5)[,2],
                    ((tabemp4.1+tabemp4.2+tabemp4.3+tabemp4.4+tabemp4.5)/5)[,2])/1000,
              cbind(((tabhype1.1+tabhype1.2+tabhype1.3+tabhype1.4+tabhype1.5)/5)[,2], #hypertension
                    ((tabhype2.1+tabhype2.2+tabhype2.3+tabhype2.4+tabhype2.5)/5)[,2],
                    ((tabhype3.1+tabhype3.2+tabhype3.3+tabhype3.4+tabhype3.5)/5)[,2],
                    ((tabhype4.1+tabhype4.2+tabhype4.3+tabhype4.4+tabhype4.5)/5)[,2])/1000,
              cbind(((tabchd1.1+tabchd1.2+tabchd1.3+tabchd1.4+tabchd1.5)/5)[,2], #chronic heart disease
                    ((tabchd2.1+tabchd2.2+tabchd2.3+tabchd2.4+tabchd2.5)/5)[,2],
                    ((tabchd3.1+tabchd3.2+tabchd3.3+tabchd3.4+tabchd3.5)/5)[,2],
                    ((tabchd4.1+tabchd4.2+tabchd4.3+tabchd4.4+tabchd4.5)/5)[,2])/1000,
              cbind(((tabang1.1+tabang1.2+tabang1.3+tabang1.4+tabang1.5)/5)[,2], #angina
                    ((tabang2.1+tabang2.2+tabang2.3+tabang2.4+tabang2.5)/5)[,2],
                    ((tabang3.1+tabang3.2+tabang3.3+tabang3.4+tabang3.5)/5)[,2],
                    ((tabang4.1+tabang4.2+tabang4.3+tabang4.4+tabang4.5)/5)[,2])/1000,
              cbind(((tabha1.1+tabha1.2+tabha1.3+tabha1.4+tabha1.5)/5)[,2], #heart attack
                    ((tabha2.1+tabha2.2+tabha2.3+tabha2.4+tabha2.5)/5)[,2],
                    ((tabha3.1+tabha3.2+tabha3.3+tabha3.4+tabha3.5)/5)[,2],
                    ((tabha4.1+tabha4.2+tabha4.3+tabha4.4+tabha4.5)/5)[,2])/1000,
              cbind(((tabhd1.1+tabhd1.2+tabhd1.3+tabhd1.4+tabhd1.5)/5)[,2], #heart disease
                    ((tabhd2.1+tabhd2.2+tabhd2.3+tabhd2.4+tabhd2.5)/5)[,2],
                    ((tabhd3.1+tabhd3.2+tabhd3.3+tabhd3.4+tabhd3.5)/5)[,2],
                    ((tabhd4.1+tabhd4.2+tabhd4.3+tabhd4.4+tabhd4.5)/5)[,2])/1000,
              cbind(((tabstroke1.1+tabstroke1.2+tabstroke1.3+tabstroke1.4+tabstroke1.5)/5)[,2], #stroke
                    ((tabstroke2.1+tabstroke2.2+tabstroke2.3+tabstroke2.4+tabstroke2.5)/5)[,2],
                    ((tabstroke3.1+tabstroke3.2+tabstroke3.3+tabstroke3.4+tabstroke3.5)/5)[,2],
                    ((tabstroke4.1+tabstroke4.2+tabstroke4.3+tabstroke4.4+tabstroke4.5)/5)[,2])/1000,
              cbind(((tabdiab1.1+tabdiab1.2+tabdiab1.3+tabdiab1.4+tabdiab1.5)/5)[,2], #diabetes
                    ((tabdiab2.1+tabdiab2.2+tabdiab2.3+tabdiab2.4+tabdiab2.5)/5)[,2],
                    ((tabdiab3.1+tabdiab3.2+tabdiab3.3+tabdiab3.4+tabdiab3.5)/5)[,2],
                    ((tabdiab4.1+tabdiab4.2+tabdiab4.3+tabdiab4.4+tabdiab4.5)/5)[,2])/1000,
              cbind(((tabbron1.1+tabbron1.2+tabbron1.3+tabbron1.4+tabbron1.5)/5)[,2], #bronchitis
                    ((tabbron2.1+tabbron2.2+tabbron2.3+tabbron2.4+tabbron2.5)/5)[,2],
                    ((tabbron3.1+tabbron3.2+tabbron3.3+tabbron3.4+tabbron3.5)/5)[,2],
                    ((tabbron4.1+tabbron4.2+tabbron4.3+tabbron4.4+tabbron4.5)/5)[,2])/1000,
              cbind(((tabkid1.1+tabkid1.2+tabkid1.3+tabkid1.4+tabkid1.5)/5)[,2], #kidney disease
                    ((tabkid2.1+tabkid2.2+tabkid2.3+tabkid2.4+tabkid2.5)/5)[,2],
                    ((tabkid3.1+tabkid3.2+tabkid3.3+tabkid3.4+tabkid3.5)/5)[,2],
                    ((tabkid4.1+tabkid4.2+tabkid4.3+tabkid4.4+tabkid4.5)/5)[,2])/1000,
              cbind(((tabliv1.1+tabliv1.2+tabliv1.3+tabliv1.4+tabliv1.5)/5)[,2], #liver disease
                    ((tabliv2.1+tabliv2.2+tabliv2.3+tabliv2.4+tabliv2.5)/5)[,2],
                    ((tabliv3.1+tabliv3.2+tabliv3.3+tabliv3.4+tabliv3.5)/5)[,2],
                    ((tabliv4.1+tabliv4.2+tabliv4.3+tabliv4.4+tabliv4.5)/5)[,2])/1000,
              cbind(((tabpc1.1+tabpc1.2+tabpc1.3+tabpc1.4+tabpc1.5)/5)[,2], #prior cancer
                    ((tabpc2.1+tabpc2.2+tabpc2.3+tabpc2.4+tabpc2.5)/5)[,2],
                    ((tabpc3.1+tabpc3.2+tabpc3.3+tabpc3.4+tabpc3.5)/5)[,2],
                    ((tabpc4.1+tabpc4.2+tabpc4.3+tabpc4.4+tabpc4.5)/5)[,2])/1000,
              cbind(((tabse1.1+tabse1.2+tabse1.3+tabse1.4+tabse1.5)/5)[,2], #spec equipment
                    ((tabse2.1+tabse2.2+tabse2.3+tabse2.4+tabse2.5)/5)[,2],
                    ((tabse3.1+tabse3.2+tabse3.3+tabse3.4+tabse3.5)/5)[,2],
                    ((tabse4.1+tabse4.2+tabse4.3+tabse4.4+tabse4.5)/5)[,2])/1000)

colnames(tab4) <- c("Both", "LCDRAT only", "LYG only", "None")
rownames(tab4) <- c("Mean LCDRAT","Mean LG, Days", 
                    "Males","Females",
                    "Whites","Blacks","Hispanics","Other",
                    "Mean Age",
                    "40-49","50-59","60-69","70-79","80-84",
                    "former","current",
                    "<10","10-<20","20-<30","30-<40","40-<50","50+",
                    "Current","0<-5","5<-10","10<-15","15<-20","20+",
                    "Mean Medical Conditions",
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

write.table(round(tab4,2),"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/quadrants.selection.csv",sep=",")