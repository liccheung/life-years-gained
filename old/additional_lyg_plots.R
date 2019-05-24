#marginal selection of USPSTF, Risk-based, and Life-years gained

rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age==85)
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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis.1 <- nhis[order(nhis$pid),]


load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age==85)
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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis.2 <- nhis[order(nhis$pid),]



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age==85)
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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis.3 <- nhis[order(nhis$pid),]



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age==85)
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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis.4 <- nhis[order(nhis$pid),]




load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age==85)
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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis.5 <- nhis[order(nhis$pid),]

nhis <- nhis[order(nhis$pid),]
nhis$lcrat <- (nhis.1$lcrat+nhis.2$lcrat+nhis.3$lcrat+nhis.4$lcrat+nhis.5$lcrat)/5
nhis$lcdeath_benefit <- (nhis.1$lcdeath_benefit+nhis.2$lcdeath_benefit+nhis.3$lcdeath_benefit+nhis.4$lcdeath_benefit+nhis.5$lcdeath_benefit)/5
nhis$lyg <- (nhis.1$lyg+nhis.2$lyg+nhis.3$lyg+nhis.4$lyg+nhis.5$lyg)/5
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcrat
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcrat.eligible <- ifelse(nhis$lcrat>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
agec2 <- svytable(~age+I(lcrat>=LCRAT.cutoff.1[2]),master)[,2]/svytable(~age,master)
agec3 <- svytable(~age+I(lyg>=lyg.cutoff.1[2]),master)[,2]/svytable(~age,master)

strattot <- function(x) {
  agesub <- subset(master,age==x&lcrat.eligible==1)
  return(svytotal(~lcdeath_benefit,agesub))
}
lcdbenefit_age_riskbased <- c(rep(0,9),sapply(seq(49,84,1),strattot))

strattot <- function(x) {
  agesub <- subset(master,age==x&lcrat.eligible==1)
  return(svytotal(~lyg,agesub))
}
lygbenefit_age_riskbased <- c(rep(0,9),sapply(seq(49,84,1),strattot))

strattot <- function(x) {
  agesub <- subset(master,age==x&lyg.eligible==1)
  return(svytotal(~lcdeath_benefit,agesub))
}
lcdbenefit_age_lgbased <- c(rep(0,9),sapply(seq(49,84,1),strattot))

strattot <- function(x) {
  agesub <- subset(master,age==x&lyg.eligible==1)
  return(svytotal(~lyg,agesub))
}
lygbenefit_age_lgbased <- c(rep(0,9),sapply(seq(49,84,1),strattot))

strattot <- function(x) {
  agesub <- subset(master,age==x)
  return(svytotal(~lcdeath_benefit,agesub))
}
lcdbenefit_age <- sapply(seq(40,84,1),strattot)

strattot <- function(x) {
  agesub <- subset(master,age==x)
  return(svytotal(~lyg,agesub))
}
lygbenefit_age <- sapply(seq(40,84,1),strattot)

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas7.jpeg",width=8,height=11.5,units='in',res=900)
par(mfrow=c(1,3))
plot(seq(40,84,1),100*agec2,type="l",lwd=2,ylim=c(0,50),yaxt="n",xlab="Age",ylab="% of ever-smokers over threshold")
lines(seq(40,84,1),100*agec3,type="l",lwd=2,col="red")
axis(side=2,at=seq(0,50,5),labels=NULL)

plot(seq(40,84,1),cumsum(lcdbenefit_age_riskbased),type="l",lwd=2,yaxt="n",xlab="Age",ylab="Cumulative preventable lung cancers")
lines(seq(40,84,1),cumsum(lcdbenefit_age_lgbased),type="l",lwd=2,col="red")
axis(side=2,at=seq(0,250000,50000),labels=NULL)

plot(seq(40,84,1),cumsum(lygbenefit_age_riskbased),type="l",lwd=2,yaxt="n",xlab="Age",ylab="Cumulative life years to be gained")
lines(seq(40,84,1),cumsum(lygbenefit_age_lgbased),type="l",lwd=2,col="red")
axis(side=2,at=seq(0,500000,100000),labels=c("0","100000","200000","300000","400000","500000"))
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas8.jpeg",width=8,height=11.5,units='in',res=900)
par(mfrow=c(1,2))
plot(seq(40,84,1),lcdbenefit_age,type="l",lwd=2,xlab="Age",ylab="Total preventable lung cancers")
plot(seq(40,84,1),lygbenefit_age,type="l",lwd=2,xlab="Age",ylab="Total life years to be gained")
axis(side=2,at=seq(0,45,5),labels=NULL)
dev.off()