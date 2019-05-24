#marginal selection of USPSTF, Risk-based, and Life-years gained

rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod1.RData")
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
nhis$lcdrat <- nhis$predict[,3]/1000
nhis$lcrat <- nhis$predict[,5]/1000
nhis$falsepos <- nhis$predict[,7]/1000
nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis.1 <- nhis[order(nhis$pid),]


load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod2.RData")
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
nhis$lcdrat <- nhis$predict[,3]/1000
nhis$lcrat <- nhis$predict[,5]/1000
nhis$falsepos <- nhis$predict[,7]/1000
nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis.2 <- nhis[order(nhis$pid),]



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod3.RData")
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
nhis$lcdrat <- nhis$predict[,3]/1000
nhis$falsepos <- nhis$predict[,7]/1000
nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis.3 <- nhis[order(nhis$pid),]



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod4.RData")
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
nhis$lcdrat <- nhis$predict[,3]/1000
nhis$falsepos <- nhis$predict[,7]/1000
nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis.4 <- nhis[order(nhis$pid),]




load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod5.RData")
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
nhis$lcdrat <- nhis$predict[,3]/1000
nhis$falsepos <- nhis$predict[,7]/1000
nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis.5 <- nhis[order(nhis$pid),]

nhis <- nhis[order(nhis$pid),]
nhis$lcrat <- (nhis.1$lcrat+nhis.2$lcrat+nhis.3$lcrat+nhis.4$lcrat+nhis.5$lcrat)/5
nhis$lcdrat <- (nhis.1$lcdrat+nhis.2$lcdrat+nhis.3$lcdrat+nhis.4$lcdrat+nhis.5$lcdrat)/5
nhis$lcdeath_benefit <- (nhis.1$lcds5+nhis.2$lcds5+nhis.3$lcds5+nhis.4$lcds5+nhis.5$lcds5)/5
nhis$lyg <- (nhis.1$lyg+nhis.2$lyg+nhis.3$lyg+nhis.4$lyg+nhis.5$lyg)/5
nhis$falsepos <- (nhis.1$falsepos+nhis.2$falsepos+nhis.3$falsepos+nhis.4$falsepos+nhis.5$falsepos)/5
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

cbind(which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0095 & 365.25*nhis$lyg<=21 & 365.25*nhis$lyg>=19),
  nhis$lcdrat[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0095 & 365.25*nhis$lyg<=21 & 365.25*nhis$lyg>=19)],
  365.25*nhis$lyg[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0095 & 365.25*nhis$lyg<=21 & 365.25*nhis$lyg>=19)],
  nhis$uspstf.eligible[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0095 & 365.25*nhis$lyg<=21 & 365.25*nhis$lyg>=19)])
cbind(which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0099 & 365.25*nhis$lyg<=10.5 & 365.25*nhis$lyg>=9.5 & nhis$uspstf.eligible==1),
      nhis$lcdrat[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0099 & 365.25*nhis$lyg<=10.5 & 365.25*nhis$lyg>=9.5 & nhis$uspstf.eligible==1)],
      365.25*nhis$lyg[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0099 & 365.25*nhis$lyg<=10.5 & 365.25*nhis$lyg>=9.5 & nhis$uspstf.eligible==1)],
      nhis$comorbidities[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0099 & 365.25*nhis$lyg<=10.5 & 365.25*nhis$lyg>=9.5 & nhis$uspstf.eligible==1)],
      nhis$race[which(nhis$lcdrat<=0.0105 & nhis$lcdrat>=0.0099 & 365.25*nhis$lyg<=10.5 & 365.25*nhis$lyg>=9.5 & nhis$uspstf.eligible==1)])
cbind(which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=11 & 365.25*nhis$lyg>=9.5),
      nhis$lcdrat[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=11 & 365.25*nhis$lyg>=9.5)],
      365.25*nhis$lyg[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=11 & 365.25*nhis$lyg>=9.5)])
cbind(which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=20.5 & 365.25*nhis$lyg>=19.5),
      nhis$lcdrat[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=20.5 & 365.25*nhis$lyg>=19.5)],
      365.25*nhis$lyg[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=20.5 & 365.25*nhis$lyg>=19.5)],
      nhis$packyears[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=20.5 & 365.25*nhis$lyg>=19.5)],
      nhis$comorbidities[which(nhis$lcdrat<=0.0205 & nhis$lcdrat>=0.0195 & 365.25*nhis$lyg<=20.5 & 365.25*nhis$lyg>=19.5)])

nhis[13542,]
nhis[39668,]
nhis[36110,]
nhis[88833,]
nhis$lcrat[c(13542,39668,36110,88833)]
1.124*nhis$lcrat[c(13542,39668,36110,88833)]
nhis$lyg[c(13542,39668,36110,88833)]/(1.124*nhis$lcrat[c(13542,39668,36110,88833)])
0.204*nhis$lcdrat[c(13542,39668,36110,88833)]
nhis$lyg[c(13542,39668,36110,88833)]/(nhis$lcdeath_benefit[c(13542,39668,36110,88833)])
nhis$falsepos[c(13542,39668,36110,88833)]