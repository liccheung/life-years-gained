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

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCRAT model (40-80, highest risk)
risk <- nhis$lcdrat
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcdrat>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

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
agec2n <- svytable(~age+I(lcdrat>=LCDRAT.cutoff.1[2]),master)[,2]
agec2d <- svytable(~age,master)
agec3n <- svytable(~age+I(lyg>=lyg.cutoff.1[2]),master)[,2]
agec3d <- svytable(~age,master)
agec4n <- svytable(~age+uspstf.eligible,master)[,2]
agec4d <- svytable(~age,master)

agec2v2 <- c(sum(agec2n[1:3])/sum(agec2d[1:3]),
             sum(agec2n[4:6])/sum(agec2d[4:6]),
             sum(agec2n[7:9])/sum(agec2d[7:9]),
             sum(agec2n[10:12])/sum(agec2d[10:12]),
             sum(agec2n[13:15])/sum(agec2d[13:15]),
             sum(agec2n[16:18])/sum(agec2d[16:18]),
             sum(agec2n[19:21])/sum(agec2d[19:21]),
             sum(agec2n[22:24])/sum(agec2d[22:24]),
             sum(agec2n[25:27])/sum(agec2d[25:27]),
             sum(agec2n[28:30])/sum(agec2d[28:30]),
             sum(agec2n[31:33])/sum(agec2d[31:33]),
             sum(agec2n[34:36])/sum(agec2d[34:36]),
             sum(agec2n[37:39])/sum(agec2d[37:39]),
             sum(agec2n[40:42])/sum(agec2d[40:42]),
             sum(agec2n[43:45])/sum(agec2d[43:45]))

agec3v2 <- c(sum(agec3n[1:3])/sum(agec3d[1:3]),
             sum(agec3n[4:6])/sum(agec3d[4:6]),
             sum(agec3n[7:9])/sum(agec3d[7:9]),
             sum(agec3n[10:12])/sum(agec3d[10:12]),
             sum(agec3n[13:15])/sum(agec3d[13:15]),
             sum(agec3n[16:18])/sum(agec3d[16:18]),
             sum(agec3n[19:21])/sum(agec3d[19:21]),
             sum(agec3n[22:24])/sum(agec3d[22:24]),
             sum(agec3n[25:27])/sum(agec3d[25:27]),
             sum(agec3n[28:30])/sum(agec3d[28:30]),
             sum(agec3n[31:33])/sum(agec3d[31:33]),
             sum(agec3n[34:36])/sum(agec3d[34:36]),
             sum(agec3n[37:39])/sum(agec3d[37:39]),
             sum(agec3n[40:42])/sum(agec3d[40:42]),
             sum(agec3n[43:45])/sum(agec3d[43:45]))

plotage <- c(40,54.9,55,57,60,63,66,69,72,75,78,80,80.1,84)
agec4v2 <- c(0,0,agec4n[16]/agec4d[16],
             sum(agec4n[17:19])/sum(agec4d[17:19]),
             sum(agec4n[20:22])/sum(agec4d[20:22]),
             sum(agec4n[23:25])/sum(agec4d[23:25]),
             sum(agec4n[26:28])/sum(agec4d[26:28]),
             sum(agec4n[29:31])/sum(agec4d[29:31]),
             sum(agec4n[32:34])/sum(agec4d[32:34]),
             sum(agec4n[35:37])/sum(agec4d[35:37]),
             sum(agec4n[38:40])/sum(agec4d[38:40]),             
             agec4n[41]/agec4d[41],0,0)

age50_70 <- subset(master,age>=50 & age<70)
age70_84 <- subset(master,age>=70)
svymean(~lcdrat.eligible,age50_70)
svymean(~lyg.eligible,age50_70)

sum(agec2n[11:20])/sum(agec2d[11:20])
sum(agec2n[21:30])/sum(agec2d[21:30])
sum(agec2n[31:45])/sum(agec2d[31:45])
sum(agec3n[11:20])/sum(agec3d[11:20])
sum(agec3n[21:30])/sum(agec3d[21:30])
sum(agec3n[31:45])/sum(agec3d[31:45])

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/lg_vs_risk.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",main="% of US ever-smokers selected for screening by age: risk-based versus life-gained-based")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=72,y=45,"Risk-based",col="red",cex=1.5)
text(x=75,y=18,"Life-gained-based",col="blue",cex=1.5)
dev.off()

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas9_v8.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",
     main="Percentage of ever-smokers selected for screening")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
lines(plotage,100*agec4v2,type="l",lwd=2)
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=73,y=45,"Risk-based",col="red",cex=1.5)
text(x=79,y=32,"Life gained \nbased",col="blue",cex=1.5)
text(x=55,y=29,"USPSTF",cex=1.5)
dev.off()

male <- subset(master,female==0)

agec2n <- svytable(~age+I(lcdrat>=LCDRAT.cutoff.1[2]),male)[,2]
agec2d <- svytable(~age,male)
agec3n <- svytable(~age+I(lyg>=lyg.cutoff.1[2]),male)[,2]
agec3d <- svytable(~age,male)
agec4n <- svytable(~age+uspstf.eligible,male)[,2]
agec4d <- svytable(~age,male)

agec2v2 <- c(sum(agec2n[1:3])/sum(agec2d[1:3]),
             sum(agec2n[4:6])/sum(agec2d[4:6]),
             sum(agec2n[7:9])/sum(agec2d[7:9]),
             sum(agec2n[10:12])/sum(agec2d[10:12]),
             sum(agec2n[13:15])/sum(agec2d[13:15]),
             sum(agec2n[16:18])/sum(agec2d[16:18]),
             sum(agec2n[19:21])/sum(agec2d[19:21]),
             sum(agec2n[22:24])/sum(agec2d[22:24]),
             sum(agec2n[25:27])/sum(agec2d[25:27]),
             sum(agec2n[28:30])/sum(agec2d[28:30]),
             sum(agec2n[31:33])/sum(agec2d[31:33]),
             sum(agec2n[34:36])/sum(agec2d[34:36]),
             sum(agec2n[37:39])/sum(agec2d[37:39]),
             sum(agec2n[40:42])/sum(agec2d[40:42]),
             sum(agec2n[43:45])/sum(agec2d[43:45]))

agec3v2 <- c(sum(agec3n[1:3])/sum(agec3d[1:3]),
             sum(agec3n[4:6])/sum(agec3d[4:6]),
             sum(agec3n[7:9])/sum(agec3d[7:9]),
             sum(agec3n[10:12])/sum(agec3d[10:12]),
             sum(agec3n[13:15])/sum(agec3d[13:15]),
             sum(agec3n[16:18])/sum(agec3d[16:18]),
             sum(agec3n[19:21])/sum(agec3d[19:21]),
             sum(agec3n[22:24])/sum(agec3d[22:24]),
             sum(agec3n[25:27])/sum(agec3d[25:27]),
             sum(agec3n[28:30])/sum(agec3d[28:30]),
             sum(agec3n[31:33])/sum(agec3d[31:33]),
             sum(agec3n[34:36])/sum(agec3d[34:36]),
             sum(agec3n[37:39])/sum(agec3d[37:39]),
             sum(agec3n[40:42])/sum(agec3d[40:42]),
             sum(agec3n[43:45])/sum(agec3d[43:45]))

plotage <- c(40,54.9,55,57,60,63,66,69,72,75,78,80,80.1,84)
agec4v2 <- c(0,0,agec4n[16]/agec4d[16],
             sum(agec4n[17:19])/sum(agec4d[17:19]),
             sum(agec4n[20:22])/sum(agec4d[20:22]),
             sum(agec4n[23:25])/sum(agec4d[23:25]),
             sum(agec4n[26:28])/sum(agec4d[26:28]),
             sum(agec4n[29:31])/sum(agec4d[29:31]),
             sum(agec4n[32:34])/sum(agec4d[32:34]),
             sum(agec4n[35:37])/sum(agec4d[35:37]),
             sum(agec4n[38:40])/sum(agec4d[38:40]),             
             agec4n[41]/agec4d[41],0,0)

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas9_v7_male.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",
     main="Percentage of male ever-smokers selected for screening")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=72,y=45,"Risk-based",col="red",cex=1.5)
text(x=75,y=18,"Life gained based",col="blue",cex=1.5)
dev.off()

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas9_v8_male.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",
     main="Percentage of male ever-smokers selected for screening")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
lines(plotage,100*agec4v2,type="l",lwd=2)
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=73,y=45,"Risk-based",col="red",cex=1.5)
text(x=79,y=32,"Life gained \nbased",col="blue",cex=1.5)
text(x=55,y=35,"USPSTF",cex=1.5)
dev.off()



female <- subset(master,female==1)

agec2n <- svytable(~age+I(lcdrat>=LCDRAT.cutoff.1[2]),female)[,2]
agec2d <- svytable(~age,female)
agec3n <- svytable(~age+I(lyg>=lyg.cutoff.1[2]),female)[,2]
agec3d <- svytable(~age,female)
agec4n <- svytable(~age+uspstf.eligible,female)[,2]
agec4d <- svytable(~age,female)

#100*sum(agec2n[11:16])/sum(agec2d[11:16])
#100*sum(agec3n[17:16])/sum(agec3d[11:16])


agec2v2 <- c(sum(agec2n[1:3])/sum(agec2d[1:3]),
             sum(agec2n[4:6])/sum(agec2d[4:6]),
             sum(agec2n[7:9])/sum(agec2d[7:9]),
             sum(agec2n[10:12])/sum(agec2d[10:12]),
             sum(agec2n[13:15])/sum(agec2d[13:15]),
             sum(agec2n[16:18])/sum(agec2d[16:18]),
             sum(agec2n[19:21])/sum(agec2d[19:21]),
             sum(agec2n[22:24])/sum(agec2d[22:24]),
             sum(agec2n[25:27])/sum(agec2d[25:27]),
             sum(agec2n[28:30])/sum(agec2d[28:30]),
             sum(agec2n[31:33])/sum(agec2d[31:33]),
             sum(agec2n[34:36])/sum(agec2d[34:36]),
             sum(agec2n[37:39])/sum(agec2d[37:39]),
             sum(agec2n[40:42])/sum(agec2d[40:42]),
             sum(agec2n[43:45])/sum(agec2d[43:45]))

agec3v2 <- c(sum(agec3n[1:3])/sum(agec3d[1:3]),
             sum(agec3n[4:6])/sum(agec3d[4:6]),
             sum(agec3n[7:9])/sum(agec3d[7:9]),
             sum(agec3n[10:12])/sum(agec3d[10:12]),
             sum(agec3n[13:15])/sum(agec3d[13:15]),
             sum(agec3n[16:18])/sum(agec3d[16:18]),
             sum(agec3n[19:21])/sum(agec3d[19:21]),
             sum(agec3n[22:24])/sum(agec3d[22:24]),
             sum(agec3n[25:27])/sum(agec3d[25:27]),
             sum(agec3n[28:30])/sum(agec3d[28:30]),
             sum(agec3n[31:33])/sum(agec3d[31:33]),
             sum(agec3n[34:36])/sum(agec3d[34:36]),
             sum(agec3n[37:39])/sum(agec3d[37:39]),
             sum(agec3n[40:42])/sum(agec3d[40:42]),
             sum(agec3n[43:45])/sum(agec3d[43:45]))

plotage <- c(40,54.9,55,57,60,63,66,69,72,75,78,80,80.1,84)
agec4v2 <- c(0,0,agec4n[16]/agec4d[16],
             sum(agec4n[17:19])/sum(agec4d[17:19]),
             sum(agec4n[20:22])/sum(agec4d[20:22]),
             sum(agec4n[23:25])/sum(agec4d[23:25]),
             sum(agec4n[26:28])/sum(agec4d[26:28]),
             sum(agec4n[29:31])/sum(agec4d[29:31]),
             sum(agec4n[32:34])/sum(agec4d[32:34]),
             sum(agec4n[35:37])/sum(agec4d[35:37]),
             sum(agec4n[38:40])/sum(agec4d[38:40]),             
             agec4n[41]/agec4d[41],0,0)

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas9_v7_female.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",
     main="Percentage of female ever-smokers selected for screening")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=72,y=45,"Risk-based",col="red",cex=1.5)
text(x=78,y=23,"Life gained based",col="blue",cex=1.5)
dev.off()

#need to add labels to whichever figure is used
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas9_v8_female.jpeg",width=8,height=11.5,units='in',res=900)
plot(seq(41.5,83.5,3),100*agec2v2,type="l",lwd=2,col="red",ylim=c(0,55),yaxt="n",xaxt="n",xlab="Age",ylab="% of ever-smokers over threshold",
     main="Percentage of female ever-smokers selected for screening")
lines(seq(41.5,83.5,3),100*agec3v2,type="l",lwd=2,col="blue")
lines(plotage,100*agec4v2,type="l",lwd=2)
axis(side=1,at=seq(40,85,5),labels=NULL)
axis(side=2,at=seq(0,55,5),labels=NULL)
text(x=73,y=45,"Risk-based",col="red",cex=1.5)
text(x=79,y=23,"Life gained \nbased",col="blue",cex=1.5)
text(x=55,y=24,"USPSTF",cex=1.5)
dev.off()