rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_1.RData")
nhis <- subset(nhis,age<=80)
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$LCDRAT <- ifelse(nhis$age<=80,nhis$lcd5,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lcd.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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


nhis$C1 <- ifelse(nhis$uspstf.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$C2 <- ifelse(nhis$uspstf.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$uspstf.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C4 <- ifelse(nhis$uspstf.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C5 <- ifelse(nhis$lyg.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C6 <- ifelse(nhis$lyg.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab1.1 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.1 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.1 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancer Deaths Saved
lcds1.1 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.1 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.1 <- svytotal(~lcd5,lyg)*(1-0.796)

#Number of Life-years Gained
lyg1.1 <- svytotal(~lyg,uspstf)
lyg2.1 <- svytotal(~lyg,lcdrat)
lyg3.1 <- svytotal(~lyg,lyg)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
lcd5.q.1 <- svyquantile(~lcd5,any,quantiles)
nhis$lcd5cat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
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
c1 <- subset(master,uspstf.eligible==1 & lcdrat.eligible==0)
c2 <- subset(master,uspstf.eligible==0 & lcdrat.eligible==1)
c3 <- subset(master,uspstf.eligible==1 & lyg.eligible==0)
c4 <- subset(master,uspstf.eligible==0 & lyg.eligible==1)
c5 <- subset(master,lyg.eligible==0 & lcdrat.eligible==1)
c6 <- subset(master,lyg.eligible==1 & lcdrat.eligible==0)

mlcd0.1 <- svymean(~lcd5,master) 
mlcd1.1 <- svymean(~lcd5,c1)
mlcd2.1 <- svymean(~lcd5,c2)
mlcd3.1 <- svymean(~lcd5,c3)
mlcd4.1 <- svymean(~lcd5,c4)
mlcd5.1 <- svymean(~lcd5,c5)
mlcd6.1 <- svymean(~lcd5,c6)
mlyg0.1 <- svymean(~lyg,master) 
mlyg1.1 <- svymean(~lyg,c1)
mlyg2.1 <- svymean(~lyg,c2)
mlyg3.1 <- svymean(~lyg,c3)
mlyg4.1 <- svymean(~lyg,c4)
mlyg5.1 <- svymean(~lyg,c5)
mlyg6.1 <- svymean(~lyg,c6)

tablcd0.1 <- svytable(~lcd5cat,master)
tablcd1.1 <- svytable(~lcd5cat+C1,master)
tablcd2.1 <- svytable(~lcd5cat+C2,master)
tablcd3.1 <- svytable(~lcd5cat+C3,master)
tablcd4.1 <- svytable(~lcd5cat+C4,master)
tablcd5.1 <- svytable(~lcd5cat+C5,master)
tablcd6.1 <- svytable(~lcd5cat+C6,master)

tablyg0.1 <- svytable(~lygcat,master)
tablyg1.1 <- svytable(~lygcat+C1,master)
tablyg2.1 <- svytable(~lygcat+C2,master)
tablyg3.1 <- svytable(~lygcat+C3,master)
tablyg4.1 <- svytable(~lygcat+C4,master)
tablyg5.1 <- svytable(~lygcat+C5,master)
tablyg6.1 <- svytable(~lygcat+C6,master)

tabg0.1 <- svytable(~female,master)
tabg1.1 <- svytable(~female+C1,master)
tabg2.1 <- svytable(~female+C2,master)
tabg3.1 <- svytable(~female+C3,master)
tabg4.1 <- svytable(~female+C4,master)
tabg5.1 <- svytable(~female+C5,master)
tabg6.1 <- svytable(~female+C6,master)

tabr0.1 <- svytable(~race,master)
tabr1.1 <- svytable(~race+C1,master)
tabr2.1 <- svytable(~race+C2,master)
tabr3.1 <- svytable(~race+C3,master)
tabr4.1 <- svytable(~race+C4,master)
tabr5.1 <- svytable(~race+C5,master)
tabr6.1 <- svytable(~race+C6,master)

meda0.1 <- svyquantile(~age,master,.1)
meda.uspstf.1 <- svyquantile(~age,uspstf,.1)
meda.lcd.1 <- svyquantile(~age,lcdrat,.1)
meda.lyg.1 <- svyquantile(~age,lyg,.1)
meda1.1 <- svyquantile(~age,c1,.1)
meda2.1 <- svyquantile(~age,c2,.1)
meda3.1 <- svyquantile(~age,c3,.1)
meda4.1 <- svyquantile(~age,c4,.1)
meda5.1 <- svyquantile(~age,c5,.1)
meda6.1 <- svyquantile(~age,c6,.1)

taba0.1 <- svytable(~age.cat,master)
taba1.1 <- svytable(~age.cat+C1,master)
taba2.1 <- svytable(~age.cat+C2,master)
taba3.1 <- svytable(~age.cat+C3,master)
taba4.1 <- svytable(~age.cat+C4,master)
taba5.1 <- svytable(~age.cat+C5,master)
taba6.1 <- svytable(~age.cat+C6,master)

tabc0.1 <- svytable(~comorbidcat,master)
tabc1.1 <- svytable(~comorbidcat+C1,master)
tabc2.1 <- svytable(~comorbidcat+C2,master)
tabc3.1 <- svytable(~comorbidcat+C3,master)
tabc4.1 <- svytable(~comorbidcat+C4,master)
tabc5.1 <- svytable(~comorbidcat+C5,master)
tabc6.1 <- svytable(~comorbidcat+C6,master)



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_2.RData")
nhis <- subset(nhis,age<=80)
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$LCDRAT <- ifelse(nhis$age<=80,nhis$lcd5,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lcd.cutoff.2 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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


nhis$C1 <- ifelse(nhis$uspstf.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$C2 <- ifelse(nhis$uspstf.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$uspstf.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C4 <- ifelse(nhis$uspstf.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C5 <- ifelse(nhis$lyg.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C6 <- ifelse(nhis$lyg.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab1.2 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.2 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.2 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancer Deaths Saved
lcds1.2 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.2 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.2 <- svytotal(~lcd5,lyg)*(1-0.796)

#Number of Life-years Gained
lyg1.2 <- svytotal(~lyg,uspstf)
lyg2.2 <- svytotal(~lyg,lcdrat)
lyg3.2 <- svytotal(~lyg,lyg)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
lcd5.q.2 <- svyquantile(~lcd5,any,quantiles)
nhis$lcd5cat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
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
c1 <- subset(master,uspstf.eligible==1 & lcdrat.eligible==0)
c2 <- subset(master,uspstf.eligible==0 & lcdrat.eligible==1)
c3 <- subset(master,uspstf.eligible==1 & lyg.eligible==0)
c4 <- subset(master,uspstf.eligible==0 & lyg.eligible==1)
c5 <- subset(master,lyg.eligible==0 & lcdrat.eligible==1)
c6 <- subset(master,lyg.eligible==1 & lcdrat.eligible==0)

mlcd0.2 <- svymean(~lcd5,master) 
mlcd1.2 <- svymean(~lcd5,c1)
mlcd2.2 <- svymean(~lcd5,c2)
mlcd3.2 <- svymean(~lcd5,c3)
mlcd4.2 <- svymean(~lcd5,c4)
mlcd5.2 <- svymean(~lcd5,c5)
mlcd6.2 <- svymean(~lcd5,c6)
mlyg0.2 <- svymean(~lyg,master) 
mlyg1.2 <- svymean(~lyg,c1)
mlyg2.2 <- svymean(~lyg,c2)
mlyg3.2 <- svymean(~lyg,c3)
mlyg4.2 <- svymean(~lyg,c4)
mlyg5.2 <- svymean(~lyg,c5)
mlyg6.2 <- svymean(~lyg,c6)

tablcd0.2 <- svytable(~lcd5cat,master)
tablcd1.2 <- svytable(~lcd5cat+C1,master)
tablcd2.2 <- svytable(~lcd5cat+C2,master)
tablcd3.2 <- svytable(~lcd5cat+C3,master)
tablcd4.2 <- svytable(~lcd5cat+C4,master)
tablcd5.2 <- svytable(~lcd5cat+C5,master)
tablcd6.2 <- svytable(~lcd5cat+C6,master)

tablyg0.2 <- svytable(~lygcat,master)
tablyg1.2 <- svytable(~lygcat+C1,master)
tablyg2.2 <- svytable(~lygcat+C2,master)
tablyg3.2 <- svytable(~lygcat+C3,master)
tablyg4.2 <- svytable(~lygcat+C4,master)
tablyg5.2 <- svytable(~lygcat+C5,master)
tablyg6.2 <- svytable(~lygcat+C6,master)

tabg0.2 <- svytable(~female,master)
tabg1.2 <- svytable(~female+C1,master)
tabg2.2 <- svytable(~female+C2,master)
tabg3.2 <- svytable(~female+C3,master)
tabg4.2 <- svytable(~female+C4,master)
tabg5.2 <- svytable(~female+C5,master)
tabg6.2 <- svytable(~female+C6,master)

tabr0.2 <- svytable(~race,master)
tabr1.2 <- svytable(~race+C1,master)
tabr2.2 <- svytable(~race+C2,master)
tabr3.2 <- svytable(~race+C3,master)
tabr4.2 <- svytable(~race+C4,master)
tabr5.2 <- svytable(~race+C5,master)
tabr6.2 <- svytable(~race+C6,master)

meda0.2 <- svyquantile(~age,master,.2)
meda.uspstf.2 <- svyquantile(~age,uspstf,.2)
meda.lcd.2 <- svyquantile(~age,lcdrat,.2)
meda.lyg.2 <- svyquantile(~age,lyg,.2)
meda1.2 <- svyquantile(~age,c1,.2)
meda2.2 <- svyquantile(~age,c2,.2)
meda3.2 <- svyquantile(~age,c3,.2)
meda4.2 <- svyquantile(~age,c4,.2)
meda5.2 <- svyquantile(~age,c5,.2)
meda6.2 <- svyquantile(~age,c6,.2)

taba0.2 <- svytable(~age.cat,master)
taba1.2 <- svytable(~age.cat+C1,master)
taba2.2 <- svytable(~age.cat+C2,master)
taba3.2 <- svytable(~age.cat+C3,master)
taba4.2 <- svytable(~age.cat+C4,master)
taba5.2 <- svytable(~age.cat+C5,master)
taba6.2 <- svytable(~age.cat+C6,master)

tabc0.2 <- svytable(~comorbidcat,master)
tabc1.2 <- svytable(~comorbidcat+C1,master)
tabc2.2 <- svytable(~comorbidcat+C2,master)
tabc3.2 <- svytable(~comorbidcat+C3,master)
tabc4.2 <- svytable(~comorbidcat+C4,master)
tabc5.2 <- svytable(~comorbidcat+C5,master)
tabc6.2 <- svytable(~comorbidcat+C6,master)



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_3.RData")
nhis <- subset(nhis,age<=80)
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$LCDRAT <- ifelse(nhis$age<=80,nhis$lcd5,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lcd.cutoff.3 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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


nhis$C1 <- ifelse(nhis$uspstf.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$C2 <- ifelse(nhis$uspstf.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$uspstf.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C4 <- ifelse(nhis$uspstf.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C5 <- ifelse(nhis$lyg.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C6 <- ifelse(nhis$lyg.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab1.3 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.3 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.3 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancer Deaths Saved
lcds1.3 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.3 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.3 <- svytotal(~lcd5,lyg)*(1-0.796)

#Number of Life-years Gained
lyg1.3 <- svytotal(~lyg,uspstf)
lyg2.3 <- svytotal(~lyg,lcdrat)
lyg3.3 <- svytotal(~lyg,lyg)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
lcd5.q.3 <- svyquantile(~lcd5,any,quantiles)
nhis$lcd5cat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
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
c1 <- subset(master,uspstf.eligible==1 & lcdrat.eligible==0)
c2 <- subset(master,uspstf.eligible==0 & lcdrat.eligible==1)
c3 <- subset(master,uspstf.eligible==1 & lyg.eligible==0)
c4 <- subset(master,uspstf.eligible==0 & lyg.eligible==1)
c5 <- subset(master,lyg.eligible==0 & lcdrat.eligible==1)
c6 <- subset(master,lyg.eligible==1 & lcdrat.eligible==0)

mlcd0.3 <- svymean(~lcd5,master) 
mlcd1.3 <- svymean(~lcd5,c1)
mlcd2.3 <- svymean(~lcd5,c2)
mlcd3.3 <- svymean(~lcd5,c3)
mlcd4.3 <- svymean(~lcd5,c4)
mlcd5.3 <- svymean(~lcd5,c5)
mlcd6.3 <- svymean(~lcd5,c6)
mlyg0.3 <- svymean(~lyg,master) 
mlyg1.3 <- svymean(~lyg,c1)
mlyg2.3 <- svymean(~lyg,c2)
mlyg3.3 <- svymean(~lyg,c3)
mlyg4.3 <- svymean(~lyg,c4)
mlyg5.3 <- svymean(~lyg,c5)
mlyg6.3 <- svymean(~lyg,c6)

tablcd0.3 <- svytable(~lcd5cat,master)
tablcd1.3 <- svytable(~lcd5cat+C1,master)
tablcd2.3 <- svytable(~lcd5cat+C2,master)
tablcd3.3 <- svytable(~lcd5cat+C3,master)
tablcd4.3 <- svytable(~lcd5cat+C4,master)
tablcd5.3 <- svytable(~lcd5cat+C5,master)
tablcd6.3 <- svytable(~lcd5cat+C6,master)

tablyg0.3 <- svytable(~lygcat,master)
tablyg1.3 <- svytable(~lygcat+C1,master)
tablyg2.3 <- svytable(~lygcat+C2,master)
tablyg3.3 <- svytable(~lygcat+C3,master)
tablyg4.3 <- svytable(~lygcat+C4,master)
tablyg5.3 <- svytable(~lygcat+C5,master)
tablyg6.3 <- svytable(~lygcat+C6,master)

tabg0.3 <- svytable(~female,master)
tabg1.3 <- svytable(~female+C1,master)
tabg2.3 <- svytable(~female+C2,master)
tabg3.3 <- svytable(~female+C3,master)
tabg4.3 <- svytable(~female+C4,master)
tabg5.3 <- svytable(~female+C5,master)
tabg6.3 <- svytable(~female+C6,master)

tabr0.3 <- svytable(~race,master)
tabr1.3 <- svytable(~race+C1,master)
tabr2.3 <- svytable(~race+C2,master)
tabr3.3 <- svytable(~race+C3,master)
tabr4.3 <- svytable(~race+C4,master)
tabr5.3 <- svytable(~race+C5,master)
tabr6.3 <- svytable(~race+C6,master)

meda0.3 <- svyquantile(~age,master,.3)
meda.uspstf.3 <- svyquantile(~age,uspstf,.3)
meda.lcd.3 <- svyquantile(~age,lcdrat,.3)
meda.lyg.3 <- svyquantile(~age,lyg,.3)
meda1.3 <- svyquantile(~age,c1,.3)
meda2.3 <- svyquantile(~age,c2,.3)
meda3.3 <- svyquantile(~age,c3,.3)
meda4.3 <- svyquantile(~age,c4,.3)
meda5.3 <- svyquantile(~age,c5,.3)
meda6.3 <- svyquantile(~age,c6,.3)

taba0.3 <- svytable(~age.cat,master)
taba1.3 <- svytable(~age.cat+C1,master)
taba2.3 <- svytable(~age.cat+C2,master)
taba3.3 <- svytable(~age.cat+C3,master)
taba4.3 <- svytable(~age.cat+C4,master)
taba5.3 <- svytable(~age.cat+C5,master)
taba6.3 <- svytable(~age.cat+C6,master)

tabc0.3 <- svytable(~comorbidcat,master)
tabc1.3 <- svytable(~comorbidcat+C1,master)
tabc2.3 <- svytable(~comorbidcat+C2,master)
tabc3.3 <- svytable(~comorbidcat+C3,master)
tabc4.3 <- svytable(~comorbidcat+C4,master)
tabc5.3 <- svytable(~comorbidcat+C5,master)
tabc6.3 <- svytable(~comorbidcat+C6,master)




load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_4.RData")
nhis <- subset(nhis,age<=80)
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$LCDRAT <- ifelse(nhis$age<=80,nhis$lcd5,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lcd.cutoff.4 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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


nhis$C1 <- ifelse(nhis$uspstf.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$C2 <- ifelse(nhis$uspstf.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$uspstf.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C4 <- ifelse(nhis$uspstf.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C5 <- ifelse(nhis$lyg.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C6 <- ifelse(nhis$lyg.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab1.4 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.4 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.4 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancer Deaths Saved
lcds1.4 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.4 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.4 <- svytotal(~lcd5,lyg)*(1-0.796)

#Number of Life-years Gained
lyg1.4 <- svytotal(~lyg,uspstf)
lyg2.4 <- svytotal(~lyg,lcdrat)
lyg3.4 <- svytotal(~lyg,lyg)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
lcd5.q.4 <- svyquantile(~lcd5,any,quantiles)
nhis$lcd5cat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
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
c1 <- subset(master,uspstf.eligible==1 & lcdrat.eligible==0)
c2 <- subset(master,uspstf.eligible==0 & lcdrat.eligible==1)
c3 <- subset(master,uspstf.eligible==1 & lyg.eligible==0)
c4 <- subset(master,uspstf.eligible==0 & lyg.eligible==1)
c5 <- subset(master,lyg.eligible==0 & lcdrat.eligible==1)
c6 <- subset(master,lyg.eligible==1 & lcdrat.eligible==0)

mlcd0.4 <- svymean(~lcd5,master) 
mlcd1.4 <- svymean(~lcd5,c1)
mlcd2.4 <- svymean(~lcd5,c2)
mlcd3.4 <- svymean(~lcd5,c3)
mlcd4.4 <- svymean(~lcd5,c4)
mlcd5.4 <- svymean(~lcd5,c5)
mlcd6.4 <- svymean(~lcd5,c6)
mlyg0.4 <- svymean(~lyg,master) 
mlyg1.4 <- svymean(~lyg,c1)
mlyg2.4 <- svymean(~lyg,c2)
mlyg3.4 <- svymean(~lyg,c3)
mlyg4.4 <- svymean(~lyg,c4)
mlyg5.4 <- svymean(~lyg,c5)
mlyg6.4 <- svymean(~lyg,c6)

tablcd0.4 <- svytable(~lcd5cat,master)
tablcd1.4 <- svytable(~lcd5cat+C1,master)
tablcd2.4 <- svytable(~lcd5cat+C2,master)
tablcd3.4 <- svytable(~lcd5cat+C3,master)
tablcd4.4 <- svytable(~lcd5cat+C4,master)
tablcd5.4 <- svytable(~lcd5cat+C5,master)
tablcd6.4 <- svytable(~lcd5cat+C6,master)

tablyg0.4 <- svytable(~lygcat,master)
tablyg1.4 <- svytable(~lygcat+C1,master)
tablyg2.4 <- svytable(~lygcat+C2,master)
tablyg3.4 <- svytable(~lygcat+C3,master)
tablyg4.4 <- svytable(~lygcat+C4,master)
tablyg5.4 <- svytable(~lygcat+C5,master)
tablyg6.4 <- svytable(~lygcat+C6,master)

tabg0.4 <- svytable(~female,master)
tabg1.4 <- svytable(~female+C1,master)
tabg2.4 <- svytable(~female+C2,master)
tabg3.4 <- svytable(~female+C3,master)
tabg4.4 <- svytable(~female+C4,master)
tabg5.4 <- svytable(~female+C5,master)
tabg6.4 <- svytable(~female+C6,master)

tabr0.4 <- svytable(~race,master)
tabr1.4 <- svytable(~race+C1,master)
tabr2.4 <- svytable(~race+C2,master)
tabr3.4 <- svytable(~race+C3,master)
tabr4.4 <- svytable(~race+C4,master)
tabr5.4 <- svytable(~race+C5,master)
tabr6.4 <- svytable(~race+C6,master)

meda0.4 <- svyquantile(~age,master,.4)
meda.uspstf.4 <- svyquantile(~age,uspstf,.4)
meda.lcd.4 <- svyquantile(~age,lcdrat,.4)
meda.lyg.4 <- svyquantile(~age,lyg,.4)
meda1.4 <- svyquantile(~age,c1,.4)
meda2.4 <- svyquantile(~age,c2,.4)
meda3.4 <- svyquantile(~age,c3,.4)
meda4.4 <- svyquantile(~age,c4,.4)
meda5.4 <- svyquantile(~age,c5,.4)
meda6.4 <- svyquantile(~age,c6,.4)

taba0.4 <- svytable(~age.cat,master)
taba1.4 <- svytable(~age.cat+C1,master)
taba2.4 <- svytable(~age.cat+C2,master)
taba3.4 <- svytable(~age.cat+C3,master)
taba4.4 <- svytable(~age.cat+C4,master)
taba5.4 <- svytable(~age.cat+C5,master)
taba6.4 <- svytable(~age.cat+C6,master)

tabc0.4 <- svytable(~comorbidcat,master)
tabc1.4 <- svytable(~comorbidcat+C1,master)
tabc2.4 <- svytable(~comorbidcat+C2,master)
tabc3.4 <- svytable(~comorbidcat+C3,master)
tabc4.4 <- svytable(~comorbidcat+C4,master)
tabc5.4 <- svytable(~comorbidcat+C5,master)
tabc6.4 <- svytable(~comorbidcat+C6,master)



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj_5.RData")
nhis <- subset(nhis,age<=80)
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
nhis$LCDRAT <- ifelse(nhis$age<=80,nhis$lcd5,0)

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lcd.cutoff.5 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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


nhis$C1 <- ifelse(nhis$uspstf.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$C2 <- ifelse(nhis$uspstf.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C3 <- ifelse(nhis$uspstf.eligible==1 & nhis$lyg.eligible==0,1,0)
nhis$C4 <- ifelse(nhis$uspstf.eligible==0 & nhis$lyg.eligible==1,1,0)
nhis$C5 <- ifelse(nhis$lyg.eligible==0 & nhis$lcdrat.eligible==1,1,0)
nhis$C6 <- ifelse(nhis$lyg.eligible==1 & nhis$lcdrat.eligible==0,1,0)
nhis$age.cat <- ifelse(nhis$age>80,4,ifelse(nhis$age>=70,3,ifelse(nhis$age>=60,2,ifelse(nhis$age>=50,1,0))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
quantiles <- c(0,.25,.5,.75,1)

uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)
any <- subset(master,uspstf.eligible|lcdrat.eligible|lyg.eligible)

#Output
#Differences in Selection
tab1.5 <- svytable(~uspstf.eligible+lcdrat.eligible,master)
tab2.5 <- svytable(~uspstf.eligible+lyg.eligible,master)
tab3.5 <- svytable(~lcdrat.eligible+lyg.eligible,master)

#Number of Lung Cancer Deaths Saved
lcds1.5 <- svytotal(~lcd5,uspstf)*(1-0.796)
lcds2.5 <- svytotal(~lcd5,lcdrat)*(1-0.796)
lcds3.5 <- svytotal(~lcd5,lyg)*(1-0.796)

#Number of Life-years Gained
lyg1.5 <- svytotal(~lyg,uspstf)
lyg2.5 <- svytotal(~lyg,lcdrat)
lyg3.5 <- svytotal(~lyg,lyg)


quantiles <- c(0,0.2,0.4,0.6,0.8,1)
lcd5.q.5 <- svyquantile(~lcd5,any,quantiles)
nhis$lcd5cat <- ifelse(nhis$lcd5<=svyquantile(~lcd5,any,quantiles)[2],1,
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
c1 <- subset(master,uspstf.eligible==1 & lcdrat.eligible==0)
c2 <- subset(master,uspstf.eligible==0 & lcdrat.eligible==1)
c3 <- subset(master,uspstf.eligible==1 & lyg.eligible==0)
c4 <- subset(master,uspstf.eligible==0 & lyg.eligible==1)
c5 <- subset(master,lyg.eligible==0 & lcdrat.eligible==1)
c6 <- subset(master,lyg.eligible==1 & lcdrat.eligible==0)


mlcd0.5 <- svymean(~lcd5,master) 
mlcd1.5 <- svymean(~lcd5,c1)
mlcd2.5 <- svymean(~lcd5,c2)
mlcd3.5 <- svymean(~lcd5,c3)
mlcd4.5 <- svymean(~lcd5,c4)
mlcd5.5 <- svymean(~lcd5,c5)
mlcd6.5 <- svymean(~lcd5,c6)
mlyg0.5 <- svymean(~lyg,master) 
mlyg1.5 <- svymean(~lyg,c1)
mlyg2.5 <- svymean(~lyg,c2)
mlyg3.5 <- svymean(~lyg,c3)
mlyg4.5 <- svymean(~lyg,c4)
mlyg5.5 <- svymean(~lyg,c5)
mlyg6.5 <- svymean(~lyg,c6)

tablcd0.5 <- svytable(~lcd5cat,master)
tablcd1.5 <- svytable(~lcd5cat+C1,master)
tablcd2.5 <- svytable(~lcd5cat+C2,master)
tablcd3.5 <- svytable(~lcd5cat+C3,master)
tablcd4.5 <- svytable(~lcd5cat+C4,master)
tablcd5.5 <- svytable(~lcd5cat+C5,master)
tablcd6.5 <- svytable(~lcd5cat+C6,master)

tablyg0.5 <- svytable(~lygcat,master)
tablyg1.5 <- svytable(~lygcat+C1,master)
tablyg2.5 <- svytable(~lygcat+C2,master)
tablyg3.5 <- svytable(~lygcat+C3,master)
tablyg4.5 <- svytable(~lygcat+C4,master)
tablyg5.5 <- svytable(~lygcat+C5,master)
tablyg6.5 <- svytable(~lygcat+C6,master)

tabg0.5 <- svytable(~female,master)
tabg1.5 <- svytable(~female+C1,master)
tabg2.5 <- svytable(~female+C2,master)
tabg3.5 <- svytable(~female+C3,master)
tabg4.5 <- svytable(~female+C4,master)
tabg5.5 <- svytable(~female+C5,master)
tabg6.5 <- svytable(~female+C6,master)

tabr0.5 <- svytable(~race,master)
tabr1.5 <- svytable(~race+C1,master)
tabr2.5 <- svytable(~race+C2,master)
tabr3.5 <- svytable(~race+C3,master)
tabr4.5 <- svytable(~race+C4,master)
tabr5.5 <- svytable(~race+C5,master)
tabr6.5 <- svytable(~race+C6,master)

meda0.5 <- svyquantile(~age,master,.5)
meda.uspstf.5 <- svyquantile(~age,uspstf,.5)
meda.lcd.5 <- svyquantile(~age,lcdrat,.5)
meda.lyg.5 <- svyquantile(~age,lyg,.5)
meda1.5 <- svyquantile(~age,c1,.5)
meda2.5 <- svyquantile(~age,c2,.5)
meda3.5 <- svyquantile(~age,c3,.5)
meda4.5 <- svyquantile(~age,c4,.5)
meda5.5 <- svyquantile(~age,c5,.5)
meda6.5 <- svyquantile(~age,c6,.5)

taba0.5 <- svytable(~age.cat,master)
taba1.5 <- svytable(~age.cat+C1,master)
taba2.5 <- svytable(~age.cat+C2,master)
taba3.5 <- svytable(~age.cat+C3,master)
taba4.5 <- svytable(~age.cat+C4,master)
taba5.5 <- svytable(~age.cat+C5,master)
taba6.5 <- svytable(~age.cat+C6,master)

tabc0.5 <- svytable(~comorbidcat,master)
tabc1.5 <- svytable(~comorbidcat+C1,master)
tabc2.5 <- svytable(~comorbidcat+C2,master)
tabc3.5 <- svytable(~comorbidcat+C3,master)
tabc4.5 <- svytable(~comorbidcat+C4,master)
tabc5.5 <- svytable(~comorbidcat+C5,master)
tabc6.5 <- svytable(~comorbidcat+C6,master)



#Average values
(lcd.cutoff.1+lcd.cutoff.2+lcd.cutoff.3+lcd.cutoff.4+lcd.cutoff.5)/5
365.25*(lyg.cutoff.1+lyg.cutoff.2+lyg.cutoff.3+lyg.cutoff.4+lyg.cutoff.5)/5
(tab1.1+tab1.2+tab1.3+tab1.4+tab1.5)/5
(tab2.1+tab2.2+tab2.3+tab2.4+tab2.5)/5
(tab3.1+tab3.2+tab3.3+tab3.4+tab3.5)/5
(lcds1.1+lcds1.2+lcds1.3+lcds1.4+lcds1.5)/5
(lcds2.1+lcds2.2+lcds2.3+lcds2.4+lcds2.5)/5
(lcds3.1+lcds3.2+lcds3.3+lcds3.4+lcds3.5)/5
(lyg1.1+lyg1.2+lyg1.3+lyg1.4+lyg1.5)/5
(lyg2.1+lyg2.2+lyg2.3+lyg2.4+lyg2.5)/5
(lyg3.1+lyg3.2+lyg3.3+lyg3.4+lyg3.5)/5

(mlcd0.1+mlcd0.2+mlcd0.3+mlcd0.4+mlcd0.5)/5
(mlcd1.1+mlcd1.2+mlcd1.3+mlcd1.4+mlcd1.5)/5
(mlcd2.1+mlcd2.2+mlcd2.3+mlcd2.4+mlcd2.5)/5
(mlcd3.1+mlcd3.2+mlcd3.3+mlcd3.4+mlcd3.5)/5
(mlcd4.1+mlcd4.2+mlcd4.3+mlcd4.4+mlcd4.5)/5
(mlcd5.1+mlcd5.2+mlcd5.3+mlcd5.4+mlcd5.5)/5
(mlcd6.1+mlcd6.2+mlcd6.3+mlcd6.4+mlcd6.5)/5
365.35*(mlyg0.1+mlyg0.2+mlyg0.3+mlyg0.4+mlyg0.5)/5
365.35*(mlyg1.1+mlyg1.2+mlyg1.3+mlyg1.4+mlyg1.5)/5
365.35*(mlyg2.1+mlyg2.2+mlyg2.3+mlyg2.4+mlyg2.5)/5
365.35*(mlyg3.1+mlyg3.2+mlyg3.3+mlyg3.4+mlyg3.5)/5
365.35*(mlyg4.1+mlyg4.2+mlyg4.3+mlyg4.4+mlyg4.5)/5
365.35*(mlyg5.1+mlyg5.2+mlyg5.3+mlyg5.4+mlyg5.5)/5
365.35*(mlyg6.1+mlyg6.2+mlyg6.3+mlyg6.4+mlyg6.5)/5
(lcd5.q.1+lcd5.q.2+lcd5.q.3+lcd5.q.4+lcd5.q.5)/5
(lyg.q.1+lyg.q.2+lyg.q.3+lyg.q.4+lyg.q.5)/5
(tablcd0.1+tablcd0.2+tablcd0.3+tablcd0.4+tablcd0.5)/5 #lcd5 overall
(tablcd1.1+tablcd1.2+tablcd1.3+tablcd1.4+tablcd1.5)/5
(tablcd2.1+tablcd2.2+tablcd2.3+tablcd2.4+tablcd2.5)/5
(tablcd3.1+tablcd3.2+tablcd3.3+tablcd3.4+tablcd3.5)/5
(tablcd4.1+tablcd4.2+tablcd4.3+tablcd4.4+tablcd4.5)/5
(tablcd5.1+tablcd5.2+tablcd5.3+tablcd5.4+tablcd5.5)/5
(tablcd6.1+tablcd6.2+tablcd6.3+tablcd6.4+tablcd6.5)/5
(tablyg0.1+tablyg0.2+tablyg0.3+tablyg0.4+tablyg0.5)/5 #lyg overall
(tablyg1.1+tablyg1.2+tablyg1.3+tablyg1.4+tablyg1.5)/5
(tablyg2.1+tablyg2.2+tablyg2.3+tablyg2.4+tablyg2.5)/5
(tablyg3.1+tablyg3.2+tablyg3.3+tablyg3.4+tablyg3.5)/5
(tablyg4.1+tablyg4.2+tablyg4.3+tablyg4.4+tablyg4.5)/5
(tablyg5.1+tablyg5.2+tablyg5.3+tablyg5.4+tablyg5.5)/5
(tablyg6.1+tablyg6.2+tablyg6.3+tablyg6.4+tablyg6.5)/5
(tabg0.1+tabg0.2+tabg0.3+tabg0.4+tabg0.5)/5 #gender
(tabg1.1+tabg1.2+tabg1.3+tabg1.4+tabg1.5)/5
(tabg2.1+tabg2.2+tabg2.3+tabg2.4+tabg2.5)/5
(tabg3.1+tabg3.2+tabg3.3+tabg3.4+tabg3.5)/5
(tabg4.1+tabg4.2+tabg4.3+tabg4.4+tabg4.5)/5
(tabg5.1+tabg5.2+tabg5.3+tabg5.4+tabg5.5)/5
(tabg6.1+tabg6.2+tabg6.3+tabg6.4+tabg6.5)/5
(tabr0.1+tabr0.2+tabr0.3+tabr0.4+tabr0.5)/5 #race
(tabr1.1+tabr1.2+tabr1.3+tabr1.4+tabr1.5)/5
(tabr2.1+tabr2.2+tabr2.3+tabr2.4+tabr2.5)/5
(tabr3.1+tabr3.2+tabr3.3+tabr3.4+tabr3.5)/5
(tabr4.1+tabr4.2+tabr4.3+tabr4.4+tabr4.5)/5
(tabr5.1+tabr5.2+tabr5.3+tabr5.4+tabr5.5)/5
(tabr6.1+tabr6.2+tabr6.3+tabr6.4+tabr6.5)/5
(meda0.1+meda0.2+meda0.3+meda0.4+meda0.5)/5 #median age
(meda.uspstf.1+meda.uspstf.2+meda.uspstf.3+meda.uspstf.4+meda.uspstf.5)/5
(meda.lcd.1+meda.lcd.2+meda.lcd.3+meda.lcd.4+meda.lcd.5)/5
(meda.lyg.1+meda.lyg.2+meda.lyg.3+meda.lyg.4+meda.lyg.5)/5
(meda1.1+meda1.2+meda1.3+meda1.4+meda1.5)/5
(meda2.1+meda2.2+meda2.3+meda2.4+meda2.5)/5
(meda3.1+meda3.2+meda3.3+meda3.4+meda3.5)/5
(meda4.1+meda4.2+meda4.3+meda4.4+meda4.5)/5
(meda5.1+meda5.2+meda5.3+meda5.4+meda5.5)/5
(meda6.1+meda6.2+meda6.3+meda6.4+meda6.5)/5
(taba0.1+taba0.2+taba0.3+taba0.4+taba0.5)/5 #age categories
(taba1.1+taba1.2+taba1.3+taba1.4+taba1.5)/5
(taba2.1+taba2.2+taba2.3+taba2.4+taba2.5)/5
(taba3.1+taba3.2+taba3.3+taba3.4+taba3.5)/5
(taba4.1+taba4.2+taba4.3+taba4.4+taba4.5)/5
(taba5.1+taba5.2+taba5.3+taba5.4+taba5.5)/5
(taba6.1+taba6.2+taba6.3+taba6.4+taba6.5)/5
(tabc0.1+tabc0.2+tabc0.3+tabc0.4+tabc0.5)/5 #comorbidities
(tabc1.1+tabc1.2+tabc1.3+tabc1.4+tabc1.5)/5
(tabc2.1+tabc2.2+tabc2.3+tabc2.4+tabc2.5)/5
(tabc3.1+tabc3.2+tabc3.3+tabc3.4+tabc3.5)/5
(tabc4.1+tabc4.2+tabc4.3+tabc4.4+tabc4.5)/5
(tabc5.1+tabc5.2+tabc5.3+tabc5.4+tabc5.5)/5
(tabc6.1+tabc6.2+tabc6.3+tabc6.4+tabc6.5)/5
