#marginal selection of USPSTF, Risk-based, and Life-years gained

rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)
library(ggplot2)

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

quantiles <- seq(0,1,.1)
lyg.q.10 <- 365.25*svyquantile(~lyg,master,quantiles)
lcdrat.q.10 <- svyquantile(~lcdrat,master,quantiles)

nhis$lcdrat.cat<-ifelse(nhis$lcdrat<lcdrat.q.10[2],1,
                       ifelse(nhis$lcdrat<lcdrat.q.10[3],2,
                              ifelse(nhis$lcdrat<lcdrat.q.10[4],3,
                                     ifelse(nhis$lcdrat<lcdrat.q.10[5],4,
                                            ifelse(nhis$lcdrat<lcdrat.q.10[6],5,
                                                   ifelse(nhis$lcdrat<lcdrat.q.10[7],6,
                                                          ifelse(nhis$lcdrat<lcdrat.q.10[8],7,
                                                                 ifelse(nhis$lcdrat<lcdrat.q.10[9],8,
                                                                        ifelse(nhis$lcdrat<lcdrat.q.10[10],9,10)))))))))

nhis$quadrant <- as.factor(ifelse(nhis$lcdrat>=0.014 & 365*nhis$lyg >= 16.2,4,
                                  ifelse(nhis$lcdrat<0.014 & 365*nhis$lyg >= 16.2,3,
                                         ifelse(nhis$lcdrat>=0.014 & 365*nhis$lyg < 16.2,1,2))))
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

100*svytable(~quadrant+uspstf.eligible,master)[,2]/svytable(~quadrant,master)
100*svytable(~lcdrat.cat+uspstf.eligible,master)[,2]/svytable(~lcdrat.cat,master)
svytable(~lcdrat.cat+I(lyg>=16.2/365.25),master)
d10_at <- subset(master,lcdrat.cat==10 & lyg>=16.2/365.25)
d10_bt <- subset(master,lcdrat.cat==10 & lyg<16.2/365.25)
svymean(~comorbidities,d10_at)
svymean(~comorbidities,d10_bt)
svytable(~comorbidities,d10_at)
svytable(~comorbidities,d10_bt)
svymean(~age,d10_at)
svymean(~age,d10_bt)
svytable(~age,d10_at)
svytable(~age,d10_bt)
svymean(~current,d10_at)
svymean(~current,d10_bt)

d9_at <- subset(master,lcdrat.cat==9 & lyg>=16.2/365.25)
d9_bt <- subset(master,lcdrat.cat==9 & lyg<16.2/365.25)
svymean(~comorbidities,d9_at)
svymean(~comorbidities,d9_bt)
svytable(~comorbidities,d9_at)
svytable(~comorbidities,d9_bt)
svymean(~age,d9_at)
svymean(~age,d9_bt)
svytable(~age,d9_at)
svytable(~age,d9_bt)
svymean(~current,d9_at)
svymean(~current,d9_bt)

d8_at <- subset(master,lcdrat.cat==8 & lyg>=16.2/365.25)
d8_bt <- subset(master,lcdrat.cat==8 & lyg<16.2/365.25)
svymean(~comorbidities,d8_at)
svymean(~comorbidities,d8_bt)
svytable(~comorbidities,d8_at)
svytable(~comorbidities,d8_bt)
svymean(~age,d8_at)
svymean(~age,d8_bt)
svytable(~age,d8_at)
svytable(~age,d8_bt)
svymean(~current,d8_at)
svymean(~current,d8_bt)

nhis <- subset(nhis, analypop==1 & age>=40 & age <= 84)

p2 <- ggplot(nhis,aes(x=100*lcdrat,y=365*lyg, color=quadrant),size=adj.wt) +
  geom_point(shape=21,alpha = 0.25) + 
  scale_x_continuous(breaks=c(0,0.61,0.99,1.4,1.81,3,5),limits=c(0,5)) +
  scale_y_continuous(breaks=c(0,10,16.2,20,30,40,50,60,70,80),limits=c(0,70)) +
  labs(x="5-year lung cancer death risk, %") + 
  labs(y="Life gained from screening, days") +
  geom_hline(yintercept=16.2,col="blue") +  
  geom_vline(xintercept=0.61,col="gray",linetype="dashed") +
  geom_vline(xintercept=0.99,col="gray",linetype="dashed") +
  geom_vline(xintercept=1.81,col="gray",linetype="dashed") +  
  geom_vline(xintercept=1.4,col="red") +
  annotate("text",x=1,y=19,label="*",size=10,fontface=2,col="gold") +
  annotate("text",x=1,y=9,label="*",size=10,fontface=2,col="gold") +  
  annotate("text",x=2,y=9,label="*",size=10,fontface=2,col="gold") +    
  annotate("text",x=2,y=19,label="*",size=10,fontface=2,col="gold") +    
  annotate("text",x=4,y=15,label="life-gained threshold",col="blue",size=7,fontface=2) +
  annotate("text",x=1.55,y=60,label="risk",col="red",size=7,fontface=2) +  
  annotate("text",x=1.75,y=57.5,label="threshold",col="red",size=7,fontface=2) +  
  annotate("text",x=0.25,y=68,label="Risk \nDeciles \n1-7",col="gray",size=6,fontface=2) +
  annotate("text",x=0.82,y=68,label="Risk \nDecile \n8",col="gray",size=6,fontface=2) +
  annotate("text",x=1.4,y=68,label="Risk \nDecile \n9",col="gray",size=6,fontface=2) +
  annotate("text",x=3,y=68,label="Risk \nDecile \n10",col="gray",size=6,fontface=2) +
  annotate("text",x=3,y=40,label="In risk-based & \nin life-gained-based selections",size=6,fontface=2) + 
  annotate("text",x=0.5,y=40,label="Not in risk-based but \nin life-gained-based selection",size=6,fontface=2) +
  annotate("text",x=3,y=2,label="In risk-based but \nnot in life-gained-based selection",size=6,fontface=2) +
  annotate("text",x=0.6,y=2,label="Not selected by risk-based \nor life-gained-based",size=6,fontface=2) +
  annotate("text",x=0.50,y=20.0,label="Ever-smoker #3",size=6,fontface=2) + 
  annotate("text",x=0.50,y=10.0,label="Ever-smoker #4",size=6,fontface=2) +   
  annotate("text",x=2.50,y=10.0,label="Ever-smoker #2",size=6,fontface=2) +   
  annotate("text",x=2.50,y=20.0,label="Ever-smoker #1",size=6,fontface=2) +  
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),legend.position="none",
        axis.title=element_text(size=18,face="bold"),
        axis.text=element_text(size=14,face="bold"))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas11.jpeg",width=11.5,height=8,units='in',res=900)
plot(p2)
dev.off()

lyg_risks <- data.frame(nhis$adj.wt,365.25*nhis$lyg,100*nhis$lcdrat)
colnames(lyg_risks) <- c("weights","days.life.gained","lung.cancer.death.risk")
write.table(lyg_risks,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/risks_lyg.csv",sep=",",row.names = FALSE)