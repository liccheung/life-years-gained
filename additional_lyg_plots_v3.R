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
lcrat.q.10 <- svyquantile(~lcrat,master,quantiles)

nhis$lcrat.cat<-ifelse(nhis$lcrat<lcrat.q.10[2],1,
                       ifelse(nhis$lcrat<lcrat.q.10[3],2,
                              ifelse(nhis$lcrat<lcrat.q.10[4],3,
                                     ifelse(nhis$lcrat<lcrat.q.10[5],4,
                                            ifelse(nhis$lcrat<lcrat.q.10[6],5,
                                                   ifelse(nhis$lcrat<lcrat.q.10[7],6,
                                                          ifelse(nhis$lcrat<lcrat.q.10[8],7,
                                                                 ifelse(nhis$lcrat<lcrat.q.10[9],8,
                                                                        ifelse(nhis$lcrat<lcrat.q.10[10],9,10)))))))))

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

svytable(~lcrat.cat+I(lyg>=16.2/365.25),master)
d10_at <- subset(master,lcrat.cat==10 & lyg>=16.2/365.25)
d10_bt <- subset(master,lcrat.cat==10 & lyg<16.2/365.25)
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

d9_at <- subset(master,lcrat.cat==9 & lyg>=16.2/365.25)
d9_bt <- subset(master,lcrat.cat==9 & lyg<16.2/365.25)
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

d8_at <- subset(master,lcrat.cat==8 & lyg>=16.2/365.25)
d8_bt <- subset(master,lcrat.cat==8 & lyg<16.2/365.25)
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


p2 <- ggplot(nhis,aes(x=100*lcrat,y=365*lyg),size=adj.wt, log10="x") +
  geom_point(shape=21,alpha = 0.25) +  
  scale_x_log10(breaks=c(0.5,1,2,3,10),limits=c(0.25,20)) +
  scale_y_continuous(limits=c(0,75)) +
  labs(x="5-year lung cancer risk, %") + 
  labs(y="Life gained from screening, days") +
  geom_hline(yintercept=16.2,col="green",linetype="dotted") +  
  geom_vline(xintercept=0.48,col="blue") +
  geom_vline(xintercept=0.69,col="blue") +
  geom_vline(xintercept=1.00,col="blue") +
  geom_vline(xintercept=1.58,col="blue") +
  geom_vline(xintercept=2.74,col="blue") +  
  geom_vline(xintercept=2.20,col="red",linetype="dotted") +    
  annotate("text",x=0.4,y=18,label="life gained based threshold=16.2 days",col="green") +
  annotate("text",x=3.5,y=0,label="risk-based threshold=2.2%",col="red") +  
  annotate("text",x=0.3,y=70,label="D1-D5") +
  annotate("text",x=0.6,y=70,label="D6") +
  annotate("text",x=0.85,y=70,label="D7") +
  annotate("text",x=1.25,y=70,label="D8") +
  annotate("text",x=2.1,y=70,label="D9") +
  annotate("text",x=5,y=70,label="D10") +
  theme(plot.title = element_text(hjust=0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas2.jpeg",width=11.5,height=8,units='in',res=900)
plot(p2)
dev.off()

lyg_risks <- data.frame(nhis$adj.wt,365.25*nhis$lyg,100*nhis$lcrat)
colnames(lyg_risks) <- c("weights","days.life.gained","lung.cancer.risk")
write.table(lyg_risks,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/risks_lyg.csv",sep=",",row.names = FALSE)